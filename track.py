import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon
import geopandas as gpd
import inspect
from netCDF4 import Dataset
import os
import pickle
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import glob
from scipy.interpolate import RegularGridInterpolator

lonmin,lonmax=140-2.5, 180+2.5
latmin,latmax=28-2.5, 40+2.5

argo_data = pd.read_parquet("Argo2014.parquet")
argo_data = argo_data.drop(columns=['Salinity_psu', 'Oxygen_flag', 'Oxygen_flag2', 'Datasets_number', 'Cycle_number', 'Float_serial_no'])

circle_enlargement_factor = 1.2  # 筛选过程中涡旋半径放大倍数
Glorys_path = '../copernicus/GLORYS'

def load_meta_data(path):
    '''
    加载meta数据，输出ACS, ACL, CS, CL四个数据集
    '''
    ACS=Dataset(os.path.join(path, 'META3.1exp_DT_allsat_Anticyclonic_short_19930101_20200307.nc'))
    ACL=Dataset(os.path.join(path, 'META3.1exp_DT_allsat_Anticyclonic_long_19930101_20200307.nc'))
    CS=Dataset(os.path.join(path, 'META3.1exp_DT_allsat_Cyclonic_short_19930101_20200307.nc'))
    CL=Dataset(os.path.join(path, 'META3.1exp_DT_allsat_Cyclonic_long_19930101_20200307.nc'))

    return ACS, ACL, CS, CL


def area_limit(DS, latmin, latmax, lonmin, lonmax):
    '''
    限制区域范围
    
    输出：涡旋序号，时间，中心点经度，中心点纬度，最值点经度，最值点纬度，边界经度，边界纬度，半径，速度边界经度，速度边界纬度
    '''
    time = DS.variables['time'][:].data
    # time = convert_date(time)
    center_lon = DS.variables['longitude'][:].data
    center_lat = DS.variables['latitude'][:].data
    max_lon = DS.variables['longitude_max'][:].data
    max_lat = DS.variables['latitude_max'][:].data
    contour_lon = DS.variables['effective_contour_longitude'][:].data
    contour_lat = DS.variables['effective_contour_latitude'][:].data
    radius = DS.variables['effective_radius'][:].data
    track = DS.variables['track'][:].data
    speed_contour_lon = DS.variables['speed_contour_longitude'][:].data
    speed_contour_lat = DS.variables['speed_contour_latitude'][:].data

    mask = (center_lat >= latmin) & (center_lat <= latmax) & \
           (center_lon >= lonmin) & (center_lon <= lonmax)
    indices = np.nonzero(mask)[0]
    ds = []
    current_track = None
    current_list = []
    for i in indices:
       if current_track is None:
            current_track = track[i]
       # 当track变化时，保存之前的列表，并新建一个列表
       if track[i] != current_track:
            ds.append(current_list)
            current_list = []
            current_track = track[i]
       current_list.append([i,time[i], center_lon[i], center_lat[i], 
              max_lon[i], max_lat[i], contour_lon[i], contour_lat[i], 
              radius[i], speed_contour_lon[i], speed_contour_lat[i]])
    if current_list:
              ds.append(current_list)
    return ds

def process_and_save(data, filename):
    processed = area_limit(data,latmin, latmax, lonmin, lonmax)
    print(f'{filename.split(".")[0].upper()} count:', len(processed))
    with open(filename, 'wb') as f:
        pickle.dump(processed, f)

def find_track(DS: list, num: int):
    '''
    寻找指定编号的涡旋轨迹
    
    轨迹中的每一个元素包含：涡旋序号，时间，中心点经度，中心点纬度，最值点经度，最值点纬度，边界经度，边界纬度，半径，速度边界经度，速度边界纬度
    '''
    for track in DS:
        for i in range(len(track)):
            if num == track[i][0]:
                return track
    raise ValueError('Track not found')

def filtered_float_data(DS: list, no: int):
    '''
    根据涡旋轨迹筛选浮标数据
    
    筛选方式：先寻找指定编号涡旋出现日期，由日期筛选浮标数据，再判断浮标是否在涡旋contour或是effective radius内。
    circle_enlargement_factor 用于扩大涡旋半径，收集effective radius理想圆外侧的部分浮标数据。
    '''
    # 获取指定编号的涡旋轨迹数据
    wanted_track = find_track(DS, no)
    t0 = np.datetime64('1950-01-01')
    num, time, center_lon, center_lat, max_lon, max_lat, contour_lon, contour_lat, radius, speed_contour_lon, speed_contour_lat = zip(*wanted_track)
    
    dates = t0 + time
    dates = pd.to_datetime(dates)
    
    argo_dates = pd.to_datetime({'year': argo_data["Year"],
                                 'month': argo_data["Month"],
                                 'day': argo_data["Day"]})
    mask = argo_dates.isin(dates)
    new_argo_data = argo_data[mask]  # 仅排除了日期不符合的数据

    needed_rows = []

    # 按 Profile_number 分组，确保每个 Profile_number 只处理一行
    for profile_number, group in new_argo_data.groupby('Profile_number'):
        # 取组中的第一行数据进行日期和地理条件的判断
        row = group.iloc[0]
        
        # 根据 Year, Month, Day 列构造日期
        row_date = pd.Timestamp(year=int(row['Year']), month=int(row['Month']), day=int(row['Day']))
        
        # 在 dates 数组中查找与 row_date 匹配的索引
        matching_indices = np.where(dates == row_date)[0]
        
        # 如果没有匹配的日期，跳过该组
        if len(matching_indices) == 0:
            continue
        
        # 取第一个匹配的索引
        i = matching_indices[0]
        
        # 基于轮廓点构造多边形
        poly_coords = list(zip(contour_lon[i], contour_lat[i]))
        poly = Polygon(poly_coords)
        pt = Point(row['Longitude'], row['Latitude'])
        inside_poly = poly.contains(pt)
        
        # 检查点是否在圆内，圆心和半径已知
        center = np.array([center_lon[i], center_lat[i]])
        point_coord = np.array([row['Longitude'], row['Latitude']])
        distance = np.linalg.norm(point_coord - center)
        inside_circle = distance <= (radius[i] / 111320) * circle_enlargement_factor
        
        # 如果点在多边形内或圆内，则将该组所有行数据加入结果中
        if inside_poly or inside_circle:
            # 添加原始索引
            group['Original_Index'] = group.index
            needed_rows.extend(group.to_dict(orient='records'))

    # 创建一个包含筛选数据的 DataFrame
    needed_data = pd.DataFrame(needed_rows)

    # 设置原始索引
    needed_data.set_index('Original_Index', inplace=True)

    # 删除显式显示的 'Original_Index' 列
    needed_data.index.name = None  # 这行删除索引名称

    return needed_data


def plot_track(DS: list, no: int):
    '''
    绘制指定编号涡旋的轨迹、相关浮标数据及其存在浮标日期的有效半径和轮廓。

    参数:
    DS (list): 涡旋轨迹数据集。
    no (int): 涡旋编号。
    '''
    wanted_track=find_track(DS, no)
    t0 = np.datetime64('1950-01-01')
    num, time, center_lon, center_lat, max_lon, max_lat, contour_lon, contour_lat, radius, speed_contour_lon, speed_contour_lat = zip(*wanted_track)
    
    argo_data_filtered = filtered_float_data(DS, no)
    needed_data = argo_data_filtered.groupby('Profile_number').apply(lambda group: group.iloc[0])
    needed_data.index.name = None
    
    dates = t0 + time
    dates = pd.to_datetime(dates)

    callers_local_vars = inspect.currentframe().f_back.f_locals.items()
    ds_names = [var_name for var_name, var_val in callers_local_vars if var_val is DS]
    if ds_names:
        ds_names = ds_names[0].upper()
    else:
        raise ValueError("UNKNOWN VARIABLE")

    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    if ds_names == 'ACS' or ds_names == 'ACL':
        colors =colors[1]
    else:
        colors =colors[0]
    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))

    fig, ax = plt.subplots(figsize=(30, 20))
    ax.set_title(f'Track {ds_names}{num[0]}', fontsize=20)
    ax.set_xlabel('Longitude', fontsize=20)
    ax.set_ylabel('Latitude', fontsize=20)
    world.plot(color='green', ax=ax)

    ax.plot(center_lon, center_lat, color=colors, label='Center Track')
    ax.plot(center_lon[0], center_lat[0], marker='o', color=colors, markersize=10)
    ax.plot(center_lon[-1], center_lat[-1], marker='x', color=colors, markersize=10)

    # 对于 needed_data 中的每个观测，按日期查找对应的涡旋轨迹点
    labeled = False
    for idx, row in needed_data.iterrows():
        obs_date = pd.Timestamp(year=int(row['Year']),
                                month=int(row['Month']),
                                day=int(row['Day']))
        matching = np.where(pd.to_datetime(dates) == obs_date)[0]
        if matching.size > 0:
            i = matching[0]
            # 以圆的形式绘制半径（将米近似转换为经纬度度数）
            if labeled == False:
                if ds_names == 'ACS' or ds_names == 'ACL':
                    circle = plt.Circle((center_lon[i], center_lat[i]),
                                radius[i] / 111320.0,
                                color='r', fill=False, linestyle='--', alpha=0.2, linewidth=1, label='Effective Radius')
                else:
                    circle = plt.Circle((center_lon[i], center_lat[i]),
                                radius[i] / 111320.0,
                                color='purple', fill=False, linestyle='--', alpha=0.2, linewidth=1, label='Effective Radius')        
            else:
                if ds_names == 'ACS' or ds_names == 'ACL':
                    circle = plt.Circle((center_lon[i], center_lat[i]),
                                    radius[i] / 111320.0,
                                    color='r', fill=False, linestyle='--', alpha=0.2, linewidth=1)
                else:
                    circle = plt.Circle((center_lon[i], center_lat[i]),
                                    radius[i] / 111320.0,
                                    color='purple', fill=False, linestyle='--', alpha=0.2, linewidth=1)
            ax.add_patch(circle)
            # 标记中心点
            # ax.plot(center_lon[i], center_lat[i], marker='o', color='black', markersize=10)
            # 绘制轮廓线
            if labeled == False:
                ax.plot(contour_lon[i], contour_lat[i], color=colors, linewidth=1, alpha=0.5, label='Effective Contour')
                labeled = True
            else:
                ax.plot(contour_lon[i], contour_lat[i], color=colors, linewidth=1, alpha=0.5)
            # 绘制日期
            dates_obs = pd.to_datetime(needed_data[['Year', 'Month', 'Day']].astype(int))
            if obs_date == dates_obs.min() or obs_date == dates_obs.max():
                ax.text(center_lon[i], center_lat[i], obs_date.strftime('%Y-%m-%d'), fontsize=20, color='black')

    # 绘制 needed_data 中的每个观测点
    ax.scatter(needed_data['Longitude'], needed_data['Latitude'], color='black', s=50, label='Floats')
    
    # 设定边界时排除META中错误的contour数据
    contour_lon_filtered = np.ma.masked_equal(contour_lon, 180.0)
    contour_lat_filtered = np.ma.masked_equal(contour_lat, 0.0)
    ax.set_xlim(np.min(contour_lon_filtered) - 0.5, np.max(contour_lon_filtered) + 0.5)
    ax.set_ylim(np.min(contour_lat_filtered) - 0.5, np.max(contour_lat_filtered) + 0.5)
    ax.set_aspect('equal')

    ax.legend(fontsize=18)

def convert_date(date):
    '''转换META数据中的日期格式'''
    t0 = np.datetime64('1950-01-01')
    date = t0 + date
    date = pd.to_datetime(date)
    return date

def plot_vertical(DS: list, no: int, show_fig: bool = False, save_fig: bool = False, color_mode: str = 'distance', variable: str = 'DO', show_colorbar: bool = False):
    '''
    根据涡旋轨迹和浮标数据，绘制浮标的垂直剖面图。

    参数:
        DS (list): 涡旋轨迹数据集。
        no (int): 涡旋编号。
        show_fig (bool): 是否显示图片，默认False。
        save_fig (bool): 是否保存图片，默认False。
        color_mode (str): 颜色模式，默认'distance'（按距离着色），可选'time'（按时间着色）。
        variable (str): 变量名称，默认'DO'（溶解氧），可选'Temp'（温度）。
        show_colorbar (bool): 是否显示颜色条，默认False。

    功能:
        对每个浮标平台，按Profile_number分组，绘制DO或Temp随深度变化的曲线。
        曲线颜色可根据浮标与涡旋中心的相对距离或采样时间（剖面号顺序）变化。
        支持图片保存与显示，可选择显示颜色条。
    '''
    original_variable_name = variable # 用于文件名
    if variable == 'DO':
        variable = 'DO_mol_kg'
    elif variable == 'Temp':
        variable = 'Temperature_degC'

    wanted_track = find_track(DS, no)
    argo_data_filtered = filtered_float_data(DS, no)

    callers_local_vars = inspect.currentframe().f_back.f_locals.items()
    ds_names = [var_name for var_name, var_val in callers_local_vars if var_val is DS]
    if ds_names:
        ds_names = ds_names[0].upper()
    else:
        raise ValueError("UNKNOWN VARIABLE: Could not automatically determine the dataset name.")
    
    if argo_data_filtered.empty:
        print(f"No Argo data found for eddy {ds_names}{no} to plot vertical profiles.")
        return

    for platform_id_val, platform_data in argo_data_filtered.groupby("Platform_number"):
        # unique_profile_numbers = platform_data['Profile_number'].unique() # 原文为 idxx_values
        # min_profile_num = min(unique_profile_numbers) # 原文为 min_idxx
        # max_profile_num = max(unique_profile_numbers) # 原文为 max_idxx
        # 使用 .agg 同时获取 min 和 max 更高效
        profile_num_agg = platform_data['Profile_number'].agg(['min', 'max'])
        min_profile_num = profile_num_agg['min']
        max_profile_num = profile_num_agg['max']
        
        fig = plt.figure(figsize=(30, 20)) # 保持大尺寸以便观察细节
        ax = plt.gca() # 获取当前axes，方便颜色条使用
        cmap = plt.cm.coolwarm

        # 收集所有剖面日期用于标题
        profile_dates_for_title = []

        for profile_num, rows in platform_data.groupby("Profile_number"):
            if rows.empty:
                continue
            
            # 提取日期，假设一个剖面内日期相同
            # current_profile_date = pd.Timestamp(year=int(rows.iloc[0]['Year']), # 原文 obs_date
            #                                     month=int(rows.iloc[0]['Month']),
            #                                     day=int(rows.iloc[0]['Day']))
            # 使用try-except确保日期转换的稳健性
            try:
                current_profile_date = pd.Timestamp(year=int(rows.iloc[0]['Year']),
                                                    month=int(rows.iloc[0]['Month']),
                                                    day=int(rows.iloc[0]['Day']))
                profile_dates_for_title.append(current_profile_date)
            except (ValueError, TypeError):
                # print(f"Skipping profile {profile_num} for platform {platform_id_val} due to invalid date.")
                continue


            color_value_normalized = 0.5 # 默认颜色值

            if color_mode == 'distance':
                # needed_data_for_dist_calc = rows.iloc[0] # 原文 needed_data
                # 确保 'Longitude' 和 'Latitude' 列存在
                if 'Longitude' not in rows.iloc[0] or 'Latitude' not in rows.iloc[0]:
                    # print(f"Skipping distance coloring for profile {profile_num} due to missing coordinates.")
                    pass # 使用默认颜色
                elif wanted_track:
                    track_dates_converted = convert_date([t[1] for t in wanted_track])
                    idx_track_list = [i for i, td in enumerate(track_dates_converted) if hasattr(td, 'date') and td.date() == current_profile_date.date()]
                    
                    if idx_track_list:
                        idx_track = idx_track_list[0]
                        center_lon, center_lat, radius = wanted_track[idx_track][2], wanted_track[idx_track][3], wanted_track[idx_track][8]
                        
                        if radius > 1e-6:
                            rel_x = (rows.iloc[0]['Longitude'] - center_lon) / (radius / 111320.0)
                            rel_y = (rows.iloc[0]['Latitude'] - center_lat) / (radius / 111320.0)
                            distance = np.sqrt(rel_x**2 + rel_y**2)
                            # normalized_distance = (distance - 0) / (1 - 0 + 1e-8) # 原文
                            color_value_normalized = 1.0 - np.clip(distance, 0.0, 1.0) # distance 已经是半径归一化的
                        # else: radius too small, use default color
                    # else: Date not found in track data, use default color
                # else: No track data, use default color
            
            elif color_mode == 'time':
                if max_profile_num > min_profile_num:
                    color_value_normalized = (profile_num - min_profile_num) / (max_profile_num - min_profile_num)
                else: # 只有一个剖面或所有剖面号相同
                    color_value_normalized = 0.0 
            
            color = cmap(np.clip(color_value_normalized, 0.0, 1.0))
            ax.plot(rows[variable], rows["Depth_m"], color=color, alpha=0.7) # label=current_profile_date.strftime("%Y-%m-%d") 移除以避免图例混乱
            
        if not profile_dates_for_title: # 如果该平台没有有效的剖面数据被绘制
            plt.close(fig)
            continue

        date_start_platform = min(profile_dates_for_title)
        date_end_platform = max(profile_dates_for_title)
        
        ax.set_ylim(-50, 2050)
        if variable == 'DO_mol_kg':
            ax.set_xlim(30,270)
        elif variable == 'Temperature_degC':
            ax.set_xlim(1, 32)
        
        ax.set_title(f"{ds_names}{no}, Platform: {int(platform_id_val)}, {date_start_platform.date()}~{date_end_platform.date()}", fontsize=20)
        ax.set_xlabel(variable, fontsize=20)
        ax.set_ylabel("Depth/m", fontsize=20)
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.invert_yaxis()
        # ax.legend() # 保持原样，通常垂直剖面图的图例会很拥挤

        if show_colorbar:
            norm_for_cbar = Normalize(vmin=0, vmax=1)
            scalar_mappable = ScalarMappable(cmap=cmap, norm=norm_for_cbar)
            scalar_mappable.set_array([]) 
            cbar = plt.colorbar(scalar_mappable, ax=ax, orientation='vertical', fraction=0.046, pad=0.04)

            if color_mode == 'distance':
                current_ticks = cbar.get_ticks()
                new_tick_labels = [f"{1.0 - t:.1f}" for t in current_ticks]
                cbar.set_ticks(current_ticks)
                cbar.set_ticklabels(new_tick_labels)
                cbar.set_label('Normalized Distance from Eddy Center (0=center, 1=edge)', fontsize=14)
            elif color_mode == 'time':
                cbar.set_label(f'Normalized Profile Sequence (Range: {min_profile_num} to {max_profile_num})', fontsize=14)

        if save_fig:
            output_dir = f"plot_vertical_{original_variable_name}"
            os.makedirs(output_dir, exist_ok=True)
            # 文件名与原函数保持一致
            file_var_short = "DO" if original_variable_name == "DO" else "T"
            plt.savefig(os.path.join(output_dir, f"{ds_names}{no}{file_var_short}{int(platform_id_val)}.png"), dpi=300, bbox_inches='tight')

        if show_fig:
            plt.show()
        plt.close(fig)
        
def plot_relative_position(DS: list, no: int, show_fig: bool = False, save_fig: bool = False, color_mode: str = 'distance', show_colorbar: bool = False):
    '''
    根据涡旋轨迹和浮标数据，绘制浮标在单位圆涡旋中的相对位置分布图。

    参数:
        DS (list): 涡旋轨迹数据集。
        no (int): 涡旋编号。
        show_fig (bool): 是否显示图片，默认False。
        save_fig (bool): 是否保存图片，默认False。
        color_mode (str): 颜色模式，默认'distance'（按距离着色），可选'time'（按剖面号顺序着色）。
        show_colorbar (bool): 是否显示颜色条，默认False。

    功能:
        对每个浮标平台，绘制其各剖面代表点相对于涡旋中心的归一化位置。
        点旁数字表示平台内剖面时序。坐标轴刻度显示真实地理坐标和相对坐标。
        观测点颜色可变。支持图片保存与显示，可选择显示颜色条。
    '''
    wanted_track = find_track(DS, no)
    argo_data_filtered = filtered_float_data(DS, no)

    callers_local_vars = inspect.currentframe().f_back.f_locals.items()
    ds_names = [var_name for var_name, var_val in callers_local_vars if var_val is DS]
    if ds_names:
        ds_names = ds_names[0].upper()
    else:
        raise ValueError("UNKNOWN VARIABLE: Could not automatically determine the dataset name.")

    if argo_data_filtered.empty:
        print(f"No Argo data found for eddy {ds_names}{no} to plot relative positions.")
        return

    for platform_id_val, platform_data in argo_data_filtered.groupby("Platform_number"):
        # 获取每个剖面的第一行数据作为代表点
        # 原文: needed_data = platform.groupby("Profile_number").apply(lambda group: group.iloc[0])
        #       needed_data.index.name = None
        # 使用 .first() 更简洁，并且可以直接用 Profile_number 作为索引或后续处理
        profile_first_rows = platform_data.groupby("Profile_number").first().reset_index() # reset_index 使 Profile_number 成为列

        if profile_first_rows.empty:
            # print(f"No profile data for platform {platform_id_val}.")
            continue

        # 准备每个剖面点对应的涡旋轨迹数据
        points_for_this_platform = []
        track_info_for_this_platform = [] # 用于计算该平台下的平均涡旋参数
        
        track_dates_converted = convert_date([t[1] for t in wanted_track]) if wanted_track else []

        for i, p_row in profile_first_rows.iterrows(): # i 将用作顺序编号
            try:
                current_date_profile = pd.Timestamp(year=int(p_row['Year']),
                                                    month=int(p_row['Month']),
                                                    day=int(p_row['Day']))
            except (ValueError, TypeError):
                # print(f"Skipping profile {p_row.get('Profile_number')} for platform {platform_id_val} due to invalid date.")
                continue

            center_lon, center_lat, radius = None, None, None
            if wanted_track:
                matches = [k for k, td in enumerate(track_dates_converted) if hasattr(td, 'date') and td.date() == current_date_profile.date()]
                if matches:
                    idx_track = matches[0]
                    center_lon = wanted_track[idx_track][2]
                    center_lat = wanted_track[idx_track][3]
                    radius = wanted_track[idx_track][8]
            
            if center_lon is not None and radius is not None and radius > 1e-6:
                 if 'Longitude' not in p_row or 'Latitude' not in p_row:
                    # print(f"Skipping point on {current_date_profile.date()} due to missing Longitude/Latitude.")
                    continue
                 rel_x = (p_row['Longitude'] - center_lon) / (radius / 111320.0)
                 rel_y = (p_row['Latitude'] - center_lat) / (radius / 111320.0)
                 points_for_this_platform.append({
                     'rel_x': rel_x, 'rel_y': rel_y, 'date': current_date_profile,
                     'profile_num_original_idx': p_row['Profile_number'], # 用于 'time' mode
                     'sequence_label': i + 1 # 平台内时序标签
                 })
                 track_info_for_this_platform.append([center_lon, center_lat, radius])
            # else: 涡旋数据不匹配或半径过小，则跳过该点
        
        if not points_for_this_platform:
            # print(f"No valid points with track data for platform {platform_id_val}.")
            continue

        fig, ax = plt.subplots(figsize=(30, 20)) # 原为 (30,20) 但相对位置图常用正方形，可考虑 (20,20) 或 (15,15)
        cmap = plt.cm.coolwarm

        # 准备 'time' 模式颜色归一化 (基于 Profile_number)
        min_prof_num_platform = min(p['profile_num_original_idx'] for p in points_for_this_platform)
        max_prof_num_platform = max(p['profile_num_original_idx'] for p in points_for_this_platform)

        for point_data in points_for_this_platform:
            rel_x = point_data['rel_x']
            rel_y = point_data['rel_y']
            
            color_value_normalized = 0.5
            if color_mode == 'distance':
                distance = np.sqrt(rel_x**2 + rel_y**2)
                color_value_normalized = 1.0 - np.clip(distance, 0.0, 1.0)
            elif color_mode == 'time':
                if max_prof_num_platform > min_prof_num_platform:
                    color_value_normalized = (point_data['profile_num_original_idx'] - min_prof_num_platform) / \
                                             (max_prof_num_platform - min_prof_num_platform)
                else:
                    color_value_normalized = 0.0
            
            color = cmap(np.clip(color_value_normalized, 0.0, 1.0))
            ax.scatter(rel_x, rel_y, color=color, s=300, zorder=5)
            ax.text(rel_x, rel_y, str(point_data['sequence_label']), weight='bold', fontsize=9, color='black', ha='center', va='center', zorder=6)
            
        ax.plot(0, 0, marker='x', color='black', markersize=16, markeredgewidth=3, label='Eddy Center (Relative)', zorder=3)
        circle = plt.Circle((0, 0), 1, color='gray', fill=False, linestyle='--', linewidth=2, label='Unit Eddy Boundary', zorder=2)
        ax.add_patch(circle)
        ax.set_aspect('equal')

        # 标题和坐标轴标签
        platform_dates = [p['date'] for p in points_for_this_platform]
        date_start_platform = min(platform_dates)
        date_end_platform = max(platform_dates)
        ax.set_title(f"{ds_names}{no}, Platform: {int(platform_id_val)}, {date_start_platform.date()}~{date_end_platform.date()}, Points: {len(points_for_this_platform)}", fontsize=20)
        ax.set_xlabel('Relative X (Eddy Radii)', fontsize=20)
        ax.set_ylabel('Relative Y (Eddy Radii)', fontsize=20)
        ax.tick_params(axis='both', which='major', labelsize=16) # 使用tick_params统一设置

        # 设置坐标轴刻度 (与 plot_relative_position_monthly 一致)
        if track_info_for_this_platform:
            mean_center_lon = np.mean([info[0] for info in track_info_for_this_platform])
            mean_center_lat = np.mean([info[1] for info in track_info_for_this_platform])
            mean_radius = np.mean([info[2] for info in track_info_for_this_platform])

            if not np.isnan(mean_center_lon) and not np.isnan(mean_center_lat) and not np.isnan(mean_radius) and mean_radius > 1e-6:
                mean_degrees = mean_radius / 111320.0
                
                tick_locs = [-1, -0.5, 0, 0.5, 1] # 使用更详细的刻度
                
                x_tick_labels = [f"{(mean_center_lon + tick_loc * mean_degrees):.2f}°\n({tick_loc:.1f})" for tick_loc in tick_locs]
                ax.set_xticks(tick_locs)
                ax.set_xticklabels(x_tick_labels)

                y_tick_labels = [f"{(mean_center_lat + tick_loc * mean_degrees):.2f}°\n({tick_loc:.1f})" for tick_loc in tick_locs]
                ax.set_yticks(tick_locs)
                ax.set_yticklabels(y_tick_labels)
        
        ax.set_xlim([-1.25, 1.25])
        ax.set_ylim([-1.25, 1.25])
        # ax.legend(fontsize=14) # 原代码legend注释掉了，保持一致

        # 添加颜色条
        if show_colorbar:
            norm_for_cbar = Normalize(vmin=0, vmax=1)
            scalar_mappable = ScalarMappable(cmap=cmap, norm=norm_for_cbar)
            scalar_mappable.set_array([])
            cbar = plt.colorbar(scalar_mappable, ax=ax, orientation='vertical', fraction=0.046, pad=0.04)

            if color_mode == 'distance':
                current_ticks = cbar.get_ticks()
                new_tick_labels = [f"{1.0 - t:.1f}" for t in current_ticks]
                cbar.set_ticks(current_ticks)
                cbar.set_ticklabels(new_tick_labels)
                cbar.set_label('Normalized Distance from Eddy Center (0=center, 1=edge)', fontsize=14)
            elif color_mode == 'time':
                cbar.set_label(f'Normalized Profile Sequence (Platform Range: {min_prof_num_platform} to {max_prof_num_platform})', fontsize=14)
        
        if save_fig:
            output_dir = "plot_relative_position"
            os.makedirs(output_dir, exist_ok=True)
            plt.savefig(os.path.join(output_dir, f"{ds_names}{no}RP{int(platform_id_val)}.png"), dpi=300, bbox_inches='tight')
        
        if show_fig:
            plt.show()
        plt.close(fig)

def plot_vertical_monthly(DS: list, no: int, month_required: list, show_fig: bool = False, save_fig: bool = False, color_mode: str = 'distance', variable: str = 'DO', show_colorbar: bool = False):
    '''
    绘制指定涡旋在指定月份内所有浮标平台的DO或Temp随深度变化的聚合垂直剖面图。

    参数:
        DS (list): 涡旋轨迹数据集。
        no (int): 涡旋编号。
        month_required (list): 需要绘制的月份列表，例如[1, 2, 3]表示1月、2月、3月。
        show_fig (bool): 是否显示图片，默认False。
        save_fig (bool): 是否保存图片，默认False。
        color_mode (str): 颜色模式，默认'distance'（按距离着色），可选'time'（按时间着色）。
        variable (str): 变量名称，默认'DO'（溶解氧），可选'Temp'（温度）。
        show_colorbar (bool): 是否显示颜色条，默认False。

    功能:
        将指定涡旋在指定月份内所有浮标平台的所有剖面数据绘制在同一张图上。
        曲线颜色可根据浮标与涡旋中心的相对距离或采样时间（全局归一化）变化。
        支持图片保存与显示，并带有颜色条。
    '''
    original_variable_name = variable # 保存原始变量名用于文件名
    if variable == 'DO':
        variable = 'DO_mol_kg'
    elif variable == 'Temp':
        variable = 'Temperature_degC'

    wanted_track = find_track(DS, no) # 获取涡旋轨迹数据
    argo_data_filtered = filtered_float_data(DS, no) # 获取与涡旋相关的浮标数据

    callers_local_vars = inspect.currentframe().f_back.f_locals.items() # 获取调用者局部变量
    ds_names = [var_name for var_name, var_val in callers_local_vars if var_val is DS]
    if ds_names:
        ds_names = ds_names[0].upper()
    else:
        raise ValueError("UNKNOWN VARIABLE")

    # 1. 预处理：收集所有符合月份要求的剖面数据
    profiles_to_plot = []
    all_profile_timestamps_for_time_mode = []
    all_profile_dates_for_title = []

    if argo_data_filtered.empty:
        print(f"No Argo data found for eddy {no} after initial filtering.")
    else:
        for _, platform_data in argo_data_filtered.groupby("Platform_number"):
            for _, profile_rows in platform_data.groupby("Profile_number"):
                if profile_rows.empty:
                    continue
                
                # 使用剖面数据的第一行来确定日期，假设一个剖面内的所有数据点日期相同
                # 确保 'Year', 'Month', 'Day' 列存在且为数值类型
                try:
                    current_date_profile = pd.Timestamp(year=int(profile_rows.iloc[0]['Year']),
                                                        month=int(profile_rows.iloc[0]['Month']),
                                                        day=int(profile_rows.iloc[0]['Day']))
                except (ValueError, TypeError) as e:
                    print(f"Skipping profile due to invalid date components: {e}. Data: {profile_rows.iloc[0].get('Year')}-{profile_rows.iloc[0].get('Month')}-{profile_rows.iloc[0].get('Day')}")
                    continue
                                
                if current_date_profile.month in month_required:
                    profiles_to_plot.append({
                        'rows': profile_rows,
                        'date': current_date_profile
                        # 'platform_id': _, # 平台ID，当前版本未使用，但保留以备将来扩展
                        # 'profile_id': __, # 剖面ID，当前版本未使用
                    })
                    all_profile_dates_for_title.append(current_date_profile)
                    if color_mode == 'time':
                        all_profile_timestamps_for_time_mode.append(current_date_profile)

    if not profiles_to_plot:
        print(f"No data found for eddy {ds_names}{no} in months {month_required}.")
        return # 如果没有数据可画，则直接返回

    # 2. 确定时间归一化的范围 (如果 color_mode == 'time') 和标题的日期范围
    min_time_for_norm = None
    max_time_for_norm = None
    if color_mode == 'time' and all_profile_timestamps_for_time_mode:
        min_time_for_norm = min(all_profile_timestamps_for_time_mode)
        max_time_for_norm = max(all_profile_timestamps_for_time_mode)
        if min_time_for_norm == max_time_for_norm and len(all_profile_timestamps_for_time_mode) > 1:
            # 如果所有时间戳相同但有多个剖面, 稍微扩展范围以避免除零，并允许颜色条显示范围
             max_time_for_norm = min_time_for_norm + pd.Timedelta(days=1) 
        elif min_time_for_norm == max_time_for_norm: # 只有一个唯一的时间戳
            # 对于单个时间点，所有剖面将获得相同的颜色（通常是色谱的起始颜色）
            pass # 不需要特殊处理，归一化会是0


    date_start_overall = min(all_profile_dates_for_title)
    date_end_overall = max(all_profile_dates_for_title)

    # 3. 开始绘图
    fig = plt.figure(figsize=(30, 20))
    cmap = plt.cm.coolwarm # 色谱

    for profile_info in profiles_to_plot:
        rows = profile_info['rows']
        current_date = profile_info['date']
        
        color_value_normalized = 0.5 # 默认归一化颜色值 (对应色谱中间)

        if color_mode == 'distance':
            needed_data = rows.iloc[0] # 假设剖面的第一行包含位置信息
            
            # 确保 'Longitude' 和 'Latitude' 列存在
            if 'Longitude' not in needed_data or 'Latitude' not in needed_data:
                print(f"Skipping profile on {current_date.date()} due to missing Longitude/Latitude. Using default color.")
            elif wanted_track: # 仅当有涡旋轨迹数据时才计算距离
                track_dates_converted = convert_date([t[1] for t in wanted_track])
                
                # 查找匹配的轨迹日期
                # 直接比较 pd.Timestamp 对象可能因时间部分不为午夜而出错，最好比较 .date()
                idx_track_list = [i for i, td in enumerate(track_dates_converted) if hasattr(td, 'date') and td.date() == current_date.date()]

                if len(idx_track_list) > 0:
                    idx_track = idx_track_list[0] # 如果一天有多个轨迹点，取第一个
                    center_lon, center_lat, radius = wanted_track[idx_track][2], wanted_track[idx_track][3], wanted_track[idx_track][8]
                    
                    # 避免半径为零或过小导致的除零错误
                    if radius <= 1e-6: # 半径过小 (例如，接近0米)
                        # print(f"Eddy radius is near zero ({radius}m) on {current_date.date()}. Distance cannot be normalized meaningfully. Using default color.")
                        distance = 0 # 可视为在中心
                        normalized_distance = 0.0
                    else:
                        rel_x = (needed_data['Longitude'] - center_lon) / (radius / 111320.0) # 使用浮点数进行除法
                        rel_y = (needed_data['Latitude'] - center_lat) / (radius / 111320.0) # 使用浮点数进行除法
                        distance = np.sqrt(rel_x**2 + rel_y**2)
                        # 归一化距离，0表示中心，1表示边缘。大于1表示在涡旋半径之外。
                        normalized_distance = distance # (distance - 0.0) / (1.0 - 0.0) # 假设distance已经是半径归一化的值
                    
                    # color = cmap(1-normalized_distance) -> 值越小（离中心近），颜色越偏向 cmap 的高端
                    color_value_normalized = 1.0 - np.clip(normalized_distance, 0.0, 1.0) # 裁剪以确保在[0,1]内，超出半径的按边缘处理
                else:
                    # print(f"Date {current_date.date()} not found in track data. Using default color for a profile.")
                    pass # color_value_normalized 保持 0.5
            else: # wanted_track 为空
                # print("No track data available for distance calculation. Using default color for all profiles.")
                pass # color_value_normalized 保持 0.5

        elif color_mode == 'time':
            if min_time_for_norm and max_time_for_norm and min_time_for_norm < max_time_for_norm:
                time_delta_total_seconds = (max_time_for_norm - min_time_for_norm).total_seconds()
                time_delta_current_seconds = (current_date - min_time_for_norm).total_seconds()
                if time_delta_total_seconds > 0: # 避免除以零
                    color_value_normalized = time_delta_current_seconds / time_delta_total_seconds
                else: # 总时间跨度为0（应该被上面的 min < max 捕获，但作为保险）
                    color_value_normalized = 0.0 # 所有剖面获得起始颜色
            elif min_time_for_norm and max_time_for_norm and min_time_for_norm == max_time_for_norm : # 所有剖面在同一时间点
                 color_value_normalized = 0.0 # 所有剖面获得起始颜色
            # else: 如果 min_time_for_norm 或 max_time_for_norm 为 None (例如，只有一个剖面且时间戳相同)，则使用默认值0.5

        color = cmap(np.clip(color_value_normalized, 0.0, 1.0)) # 再次确保裁剪， cmap 输入应在 [0,1]
        
        plt.plot(rows[variable], rows["Depth_m"], color=color, alpha=0.7)

    # 4. 设置图像属性和颜色条
    plt.ylim(-50, 2050) # 深度范围
    if variable == 'DO_mol_kg':
        plt.xlim(30, 270) # DO范围
    elif variable == 'Temperature_degC':
        plt.xlim(1, 32) # 温度范围
    
    month_required_str = ", ".join(map(str, month_required)) # 将月份列表转换为字符串
    title_str = f"{ds_names}{no}, Months: {month_required_str}"
    if all_profile_dates_for_title: # 仅当有实际绘制的日期时添加日期范围
        title_str += f", Data: {date_start_overall.date()}~{date_end_overall.date()}"
    plt.title(title_str, fontsize=20)
    plt.xlabel(variable, fontsize=20)
    plt.ylabel("Depth/m", fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=16)
    plt.gca().invert_yaxis() # 反转Y轴，深度向下增加

    # 添加颜色条
    if show_colorbar:
        norm_for_cbar = Normalize(vmin=0, vmax=1) 
        scalar_mappable = ScalarMappable(cmap=cmap, norm=norm_for_cbar)
        scalar_mappable.set_array([]) 
        cbar = plt.colorbar(scalar_mappable, ax=plt.gca())

        if color_mode == 'distance':
            # 颜色条刻度默认对应 (1 - Normalized Distance)。
            # 本处更新刻度标签，使其直接显示 Normalized Distance (0.0=中心, 1.0=边缘)。
            
            current_ticks = cbar.get_ticks() # 获取当前刻度位置 (其值代表 1-dist)
            
            new_tick_labels = []
            for t_val_one_minus_dist in current_ticks:
                dist_val = 1.0 - t_val_one_minus_dist # 计算对应的 normalized_distance
                label = f"{dist_val:.1f}" # 格式化为一位小数 (例如 "0.0", "0.5", "1.0")
                new_tick_labels.append(label)
                
            cbar.set_ticks(current_ticks) # 设置刻度位置
            cbar.set_ticklabels(new_tick_labels) # 设置新的刻度标签
            cbar.set_label('Normalized Distance from Eddy Center (0=center, 1=edge)', fontsize=14)
            # 此时，颜色条蓝色端 (dist=1) 标签为 "1.0", 红色端 (dist=0) 标签为 "0.0"
        elif color_mode == 'time':
            if min_time_for_norm and max_time_for_norm:
                cbar.set_label(f'Color mapped by Normalized Time ({min_time_for_norm.strftime("%Y-%m-%d")} to {max_time_for_norm.strftime("%Y-%m-%d")})', fontsize=14)
            else:
                cbar.set_label('Color mapped by Normalized Time (single or undefined time range)', fontsize=14)


    # 5. 保存和显示图片
    if save_fig:
        output_dir = f"plot_vertical_{original_variable_name}_monthly_aggregated" # 修改目录名以区分
        os.makedirs(output_dir, exist_ok=True)
        
        # 生成更合适的文件名
        month_suffix = "_".join(map(str, month_required))
        base_filename = f"{ds_names}{no}_{original_variable_name}_months_{month_suffix}_aggregated.png"
        
        plt.savefig(os.path.join(output_dir, base_filename), dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {os.path.join(output_dir, base_filename)}")

    if show_fig:
        plt.show()
    
    plt.close(fig) # 关闭图像，释放内存
    
def plot_relative_position_monthly(DS: list, no: int, month_required: list, show_fig: bool = False, save_fig: bool = False, color_mode: str = 'distance', show_colorbar: bool = False):
    '''
    根据涡旋轨迹和浮标数据，在指定月份绘制所有浮标在单位圆涡旋中的聚合相对位置分布图。

    参数:
        DS (list): 涡旋轨迹数据集。
        no (int): 涡旋编号。
        month_required (list): 需要绘制的月份列表，例如[7, 8]表示7月、8月。
        show_fig (bool): 是否显示图片，默认False。
        save_fig (bool): 是否保存图片，默认False。
        color_mode (str): 颜色模式，默认'distance'（按距离着色），可选'time'（按时间着色）。
        show_colorbar (bool): 是否显示颜色条，默认False。

    功能:
        对所有浮标平台，筛选指定月份内的剖面数据，计算每个剖面代表点相对于涡旋中心的归一化位置。
        所有点绘制在同一张图上。点旁边的数字根据日期在所选月份范围内的连续天数编号。
        观测点颜色可根据与涡旋中心的相对距离或全局采样时间变化。
        绘制单位圆表示涡旋边界。坐标轴刻度反映总体平均的真实地理坐标。
    '''
    wanted_track = find_track(DS, no)
    argo_data_filtered = filtered_float_data(DS, no)

    callers_local_vars = inspect.currentframe().f_back.f_locals.items()
    ds_names = [var_name for var_name, var_val in callers_local_vars if var_val is DS]
    if ds_names:
        ds_names = ds_names[0].upper()
    else:
        raise ValueError("UNKNOWN VARIABLE")

    # 1. 预处理：收集所有符合月份要求的剖面数据点及其相关信息
    points_to_process = [] # 存储 (日期, 剖面首行数据)
    if not argo_data_filtered.empty:
        for _, platform_data in argo_data_filtered.groupby("Platform_number"):
            # 获取每个剖面的第一行数据作为代表点
            profile_first_rows = platform_data.groupby("Profile_number").first().reset_index()
            for _, p_row in profile_first_rows.iterrows():
                try:
                    current_date_profile = pd.Timestamp(year=int(p_row['Year']),
                                                        month=int(p_row['Month']),
                                                        day=int(p_row['Day']))
                except (ValueError, TypeError) as e:
                    # print(f"Skipping profile point due to invalid date: {e}")
                    continue
                
                if current_date_profile.month in month_required:
                    points_to_process.append({'date': current_date_profile, 'data_row': p_row})
    
    if not points_to_process:
        print(f"No data found for eddy {ds_names}{no} in months {month_required}.")
        return

    # 确定日期编号的参考起始日期
    # 参考日是所选月份中最早月份的第一天，年份取实际绘制数据中的最早年份
    min_plot_date_overall = min(p['date'] for p in points_to_process)
    reference_start_date_for_labels = pd.Timestamp(year=min_plot_date_overall.year, 
                                                   month=min(month_required), 
                                                   day=1)

    points_to_plot = []
    all_track_info_for_overall_mean = [] # 用于计算平均涡旋参数
    all_profile_dates_for_title = []
    all_profile_timestamps_for_time_mode = []
    
    track_dates_converted = convert_date([t[1] for t in wanted_track]) if wanted_track else []

    for point_info in points_to_process:
        current_date = point_info['date']
        p_row = point_info['data_row']
        
        day_label = (current_date - reference_start_date_for_labels).days + 1
        
        # 匹配涡旋轨迹数据
        center_lon, center_lat, radius = None, None, None
        if wanted_track:
            matches = [i for i, td in enumerate(track_dates_converted) if hasattr(td, 'date') and td.date() == current_date.date()]
            if matches:
                idx_track = matches[0]
                center_lon = wanted_track[idx_track][2]
                center_lat = wanted_track[idx_track][3]
                radius = wanted_track[idx_track][8]

        if center_lon is not None and radius is not None and radius > 1e-6:
            # 确保 'Longitude' 和 'Latitude' 列存在于 p_row
            if 'Longitude' not in p_row or 'Latitude' not in p_row:
                # print(f"Skipping point on {current_date.date()} due to missing Longitude/Latitude.")
                continue

            rel_x = (p_row['Longitude'] - center_lon) / (radius / 111320.0)
            rel_y = (p_row['Latitude'] - center_lat) / (radius / 111320.0)
            
            points_to_plot.append({
                'rel_x': rel_x, 'rel_y': rel_y, 
                'date': current_date, 'day_label': day_label
            })
            all_track_info_for_overall_mean.append([center_lon, center_lat, radius])
            all_profile_dates_for_title.append(current_date)
            if color_mode == 'time':
                all_profile_timestamps_for_time_mode.append(current_date)
        # else: 涡旋数据不匹配或半径过小，则跳过该点

    if not points_to_plot:
        print(f"No valid points with track data found for eddy {ds_names}{no} in months {month_required}.")
        return

    # 2. 确定时间归一化的范围 (如果 color_mode == 'time') 和标题的日期范围
    min_time_for_norm = None
    max_time_for_norm = None
    if color_mode == 'time' and all_profile_timestamps_for_time_mode:
        min_time_for_norm = min(all_profile_timestamps_for_time_mode)
        max_time_for_norm = max(all_profile_timestamps_for_time_mode)
        if min_time_for_norm == max_time_for_norm and len(all_profile_timestamps_for_time_mode) > 1:
             max_time_for_norm = min_time_for_norm + pd.Timedelta(days=1) 

    date_start_overall = min(all_profile_dates_for_title)
    date_end_overall = max(all_profile_dates_for_title)

    # 3. 开始绘图
    fig, ax = plt.subplots(figsize=(30, 20)) # 保持与原函数一致的大尺寸
    cmap = plt.cm.coolwarm

    for point in points_to_plot:
        rel_x = point['rel_x']
        rel_y = point['rel_y']
        current_date = point['date']
        day_label = point['day_label']
        
        color_value_normalized = 0.5 

        if color_mode == 'distance':
            distance_from_center = np.sqrt(rel_x**2 + rel_y**2)
            # normalized_distance 的范围是 [0, ~1], 1 表示在涡旋边缘
            # cmap(1-normalized_distance) 使中心点 (dist=0) 颜色值高 (如红色), 边缘点 (dist=1) 颜色值低 (如蓝色)
            color_value_normalized = 1.0 - np.clip(distance_from_center, 0.0, 1.0)
        elif color_mode == 'time':
            if min_time_for_norm and max_time_for_norm and min_time_for_norm < max_time_for_norm:
                time_delta_total = (max_time_for_norm - min_time_for_norm).total_seconds()
                time_delta_current = (current_date - min_time_for_norm).total_seconds()
                if time_delta_total > 0:
                    color_value_normalized = time_delta_current / time_delta_total
                else: 
                    color_value_normalized = 0.0 
            elif min_time_for_norm and max_time_for_norm and min_time_for_norm == max_time_for_norm :
                 color_value_normalized = 0.0
        
        color = cmap(np.clip(color_value_normalized, 0.0, 1.0))
        
        ax.scatter(rel_x, rel_y, color=color, s=300, zorder=5) # zorder确保点在圆和叉上方
        ax.text(rel_x, rel_y, str(day_label), weight='bold', fontsize=9, color='black', ha='center', va='center', zorder=6)

    # 绘制涡旋中心点和单位圆 (只绘制一次)
    ax.plot(0, 0, marker='x', color='black', markersize=16, markeredgewidth=3, label='Eddy Center (Relative)', zorder=3)
    circle = plt.Circle((0, 0), 1, color='gray', fill=False, linestyle='--', linewidth=2, label='Unit Eddy Boundary', zorder=2)
    ax.add_patch(circle)
    ax.set_aspect('equal')

    # 4. 设置标题和坐标轴
    month_required_str = ", ".join(map(str, month_required))
    title_str = f"{ds_names}{no}, Months: {month_required_str}, Relative Positions"
    if all_profile_dates_for_title:
        title_str += f"\nData: {date_start_overall.date()}~{date_end_overall.date()}, Total Points: {len(points_to_plot)}"
    ax.set_title(title_str, fontsize=20)
    
    ax.set_xlabel('Relative X (Eddy Radii)', fontsize=20) # 标签改为相对单位
    ax.set_ylabel('Relative Y (Eddy Radii)', fontsize=20) # 标签改为相对单位
    plt.tick_params(axis='both', which='major', labelsize=16) # 调整刻度字体大小

    # 设置坐标轴刻度以反映平均真实地理坐标
    if all_track_info_for_overall_mean:
        mean_center_lon = np.mean([info[0] for info in all_track_info_for_overall_mean])
        mean_center_lat = np.mean([info[1] for info in all_track_info_for_overall_mean])
        mean_radius = np.mean([info[2] for info in all_track_info_for_overall_mean])

        if not np.isnan(mean_center_lon) and not np.isnan(mean_center_lat) and not np.isnan(mean_radius) and mean_radius > 1e-6:
            mean_degrees = mean_radius / 111320.0 # 1度约111.32公里
            
            x_tick_locs = [-1, -0.5, 0, 0.5, 1] # 相对坐标刻度位置
            x_tick_labels = [f"{(mean_center_lon + tick_loc * mean_degrees):.2f}°\n({tick_loc})" for tick_loc in x_tick_locs]
            ax.set_xticks(x_tick_locs)
            ax.set_xticklabels(x_tick_labels) # 默认字体大小已通过tick_params设置

            y_tick_locs = [-1, -0.5, 0, 0.5, 1]
            y_tick_labels = [f"{(mean_center_lat + tick_loc * mean_degrees):.2f}°\n({tick_loc})" for tick_loc in y_tick_locs]
            ax.set_yticks(y_tick_locs)
            ax.set_yticklabels(y_tick_labels)
            
    ax.set_xlim([-1.25, 1.25])
    ax.set_ylim([-1.25, 1.25]) # 修正：应为 set_ylim

    # ax.legend(fontsize=14) # 根据需要取消注释图例

    # 5. 添加颜色条
    if show_colorbar:
        norm_for_cbar = Normalize(vmin=0, vmax=1)
        scalar_mappable = ScalarMappable(cmap=cmap, norm=norm_for_cbar)
        scalar_mappable.set_array([])
        cbar = plt.colorbar(scalar_mappable, ax=ax, orientation='vertical', fraction=0.046, pad=0.04) # 调整颜色条大小和位置

        if color_mode == 'distance':
            current_ticks = cbar.get_ticks()
            new_tick_labels = [f"{1.0 - t:.1f}" for t in current_ticks]
            # 确保0和1的标签是 "0.0" 和 "1.0"
            new_tick_labels = [ "0.0" if lbl == "-0.0" else ("1.0" if lbl == "1.0" and (1.0 - current_ticks[i]) < 0.01 else lbl) for i, lbl in enumerate(new_tick_labels)]
            # 更简洁的格式化，确保0和1是一位小数
            new_tick_labels = []
            for t_val_one_minus_dist in current_ticks:
                dist_val = 1.0 - t_val_one_minus_dist
                label = f"{dist_val:.1f}"
                new_tick_labels.append(label)

            cbar.set_ticks(current_ticks)
            cbar.set_ticklabels(new_tick_labels)
            cbar.set_label('Normalized Distance from Eddy Center (0=center, 1=edge)', fontsize=14)
        elif color_mode == 'time':
            if min_time_for_norm and max_time_for_norm:
                cbar.set_label(f'Normalized Time ({min_time_for_norm.strftime("%Y-%m-%d")} to {max_time_for_norm.strftime("%Y-%m-%d")})', fontsize=14)
            else:
                cbar.set_label('Normalized Time', fontsize=14)

    # 6. 保存和显示图片
    if save_fig:
        output_dir = "plot_relative_position_monthly_aggregated" # 修改目录名
        os.makedirs(output_dir, exist_ok=True)
        month_suffix = "_".join(map(str, month_required))
        base_filename = f"{ds_names}{no}_RP_months_{month_suffix}_aggregated.png" # 修改文件名
        plt.savefig(os.path.join(output_dir, base_filename), dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {os.path.join(output_dir, base_filename)}")

    if show_fig:
        plt.show()
    
    plt.close(fig)

def get_glorys_filepath(date)-> str:
    '''
    根据给定日期返回对应的GLORYS NetCDF文件路径。

    构造目标日期的文件夹路径和文件名模式，查找匹配的文件。如果找到且仅有一个匹配文件，则返回其完整路径；
    如果没有找到或找到多个匹配文件，则抛出异常。

    参数:
        date: 由convert_date得出的需要查找的日期。

    返回:
        str: 匹配的GLORYS NetCDF文件的完整路径。

    异常:
        RuntimeError: 如果找到多个匹配文件。
        FileNotFoundError: 如果没有找到匹配文件。
    '''
    nc_path = os.path.join(Glorys_path, date.strftime("%Y"), date.strftime("%m"), 'mercatorglorys12v1_gl12_mean_')
    nc_path += date.strftime("%Y%m%d")
    matches = glob.glob(nc_path + '_R*')
    if matches:
        if len(matches) > 1:
            raise RuntimeError(f"Multiple matching files found for {date.strftime('%Y-%m-%d')}: {matches}")
        nc_path = matches[0]
        #print("Matching file:", nc_path)
        return nc_path
    else:
        raise FileNotFoundError(f"No matching file found for {date.strftime('%Y-%m-%d')}.")
        
def print_glorys_variable(nc_path: str):
    '''
    打印指定路径对应NetCDF文件中的所有变量名称、标准名称、维度和形状信息。
    '''
    nc_file = Dataset(nc_path, 'r')
    
    with Dataset(nc_path, 'r') as nc_file:
        for var in nc_file.variables:
            print(var, nc_file.variables[var].standard_name)
            print(nc_file.variables[var].dimensions, nc_file.variables[var].shape)

def find_track_glorys_filepath(DS: list, no: int) -> dict:
    '''
    根据涡旋编号在GLORYS数据集中查找对应的轨迹文件路径。

    参数:
        DS (list): 涡旋轨迹数据集。
        no (int): 涡旋编号。

    返回:
        dict: 涡旋轨迹数据文件路径的字典，格式为 {date: glorys_filepath}。
              如果未找到对应的轨迹数据或文件路径，则返回空字典。
    '''
    wanted_track = find_track(DS, no)
    if not wanted_track:
        print(f"未找到涡旋 {no} 的轨迹数据。")
        return {}
    
    glorys_filepaths_dict = {}
    for track_point in wanted_track:
        try:
            date = convert_date(track_point[1])  
            glorys_filepath = get_glorys_filepath(date)
            glorys_filepaths_dict[date] = glorys_filepath          
        except (RuntimeError, FileNotFoundError) as e:
            print(f"为涡旋 {no} 在日期 {track_point[1]} (转换后: {date if 'date' in locals() else '未知'}) 查找 GLORYS 文件时出错: {e}")
        except IndexError:
            print(f"涡旋 {no} 的轨迹点数据格式不正确: {track_point}")
        except Exception as e: # 捕获其他可能的 convert_date 或 get_glorys_filepath 异常
            print(f"处理涡旋 {no} 在日期 {track_point[1]} 时发生未知错误: {e}")

    if not glorys_filepaths_dict:
        print(f"未找到涡旋 {no} 的 GLORYS 文件。")
        return {}
        
    return glorys_filepaths_dict

def plot_track_area_surface_glorys(DS: list, no: int, needed_idx: int, variable: str = 'vorticity',
                                   show_fig: bool = False, save_fig: bool = False, deep_argo: bool = False,
                                   k: float = None, b: float = None):
    '''
    绘制指定涡旋在特定时刻的表层物理场快照及相关的Argo浮标数据。

    该函数会生成一张综合图，展示单个涡旋在某一天的详细情况。图中包括了
    GLORYS数据的表层物理场作为背景，涡旋的完整轨迹、当前位置、轮廓和
    半径，以及在该区域内符合条件的Argo浮标位置。

    参数:
        DS (list): 包含所有涡旋轨迹信息的数据集。
        no (int): 需要绘制的涡旋的唯一编号。
        needed_idx (int): 涡旋轨迹的时间点索引，用于确定绘图的具体日期。
        variable (str): 作为背景场绘制的GLORYS物理变量。默认为 'vorticity'。
        show_fig (bool): 是否在运行时显示生成的图像。默认为 False。
        save_fig (bool): 是否将生成的图像保存为文件。默认为 False。
        deep_argo (bool): 是否使用深层Argo数据模式。若为 True，则筛选700m深度的Argo数据并按溶解氧着色；若为 False，则使用表层数据。默认为 False。
        k (float, optional): 直线方程 y = kx + b 中的斜率。默认为 None。
        b (float, optional): 直线方程 y = kx + b 中的截距。默认为 None。
    '''
    wanted_track = find_track(DS, no)
    num, time, center_lon, center_lat, max_lon, max_lat, contour_lon, contour_lat, radius, speed_contour_lon, speed_contour_lat = zip(*wanted_track)
    dates = convert_date(time) if time else None

    # 获取Argo浮标数据
    argo_data_filtered = filtered_float_data(DS, no)
    argo_data_filtered = argo_data_filtered[pd.to_datetime(argo_data_filtered[['Year', 'Month', 'Day']])==dates[needed_idx]]    # 日期筛选
    if deep_argo:
        filtered_by_depth = argo_data_filtered[argo_data_filtered['Depth_m'] >= 500].copy()
        if filtered_by_depth.empty:
            print("Warning: No data found with Depth_m >= 500.")
            needed_data = pd.DataFrame(columns=argo_data_filtered.columns) # 返回一个空的DataFrame
        else:
            idx_max_do = filtered_by_depth.groupby('Profile_number')['DO_mol_kg'].idxmax()
            needed_data = filtered_by_depth.loc[idx_max_do]
            needed_data.index.name = None

    else:
        needed_data = argo_data_filtered.groupby('Profile_number').apply(lambda group: group.iloc[0])
        needed_data.index.name = None

    # 获取区域边界
    contour_lon_filtered = np.ma.masked_equal(contour_lon, 180.0)
    contour_lat_filtered = np.ma.masked_equal(contour_lat, 0.0)

    glorys_lon_min = np.min(contour_lon_filtered) - 0.5
    glorys_lon_max = np.max(contour_lon_filtered) + 0.5
    glorys_lat_min = np.min(contour_lat_filtered) - 0.5
    glorys_lat_max = np.max(contour_lat_filtered) + 0.5

    #获取背景场数据
    if variable == 'vorticity':
        glorys_lon_filtered, glorys_lat_filtered, glorys_depth_filtered, glorys_variables_filtered = get_track_area_glorys(DS, no, needed_idx, variables=['u', 'v'], surface=True)
        zeta, f = calculate_vorticity(glorys_lon_filtered, glorys_lat_filtered, glorys_variables_filtered['u'], glorys_variables_filtered['v'])
        glorys_variable_filtered = zeta/f
    else:
        glorys_lon_filtered, glorys_lat_filtered, glorys_depth_filtered, glorys_variables_filtered = get_track_area_glorys(DS, no, needed_idx, variables=[variable], surface=True)
        glorys_variable_filtered = glorys_variables_filtered[variable]

    callers_local_vars = inspect.currentframe().f_back.f_locals.items()
    ds_names = [var_name for var_name, var_val in callers_local_vars if var_val is DS]
    if ds_names:
        ds_names = ds_names[0].upper()
    else:
        raise ValueError("UNKNOWN VARIABLE")

    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    if ds_names == 'ACS' or ds_names == 'ACL':
        colors =colors[1]
    else:
        colors =colors[0]
    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))

    fig, ax = plt.subplots(figsize=(30, 25))
    ax.set_title(f'Track {ds_names}{num[0]}, {dates[needed_idx].strftime('%Y-%m-%d')}', fontsize=20)
    ax.set_xlabel('Longitude', fontsize=20)
    ax.set_ylabel('Latitude', fontsize=20)
    world.plot(color='green', ax=ax)

    # 绘制涡旋轨迹
    ax.plot(center_lon, center_lat, color=colors, label='Center Track')
    ax.plot(center_lon[0], center_lat[0], marker='o', color=colors, markersize=10)
    ax.plot(center_lon[-1], center_lat[-1], marker='x', color=colors, markersize=10)

    # 绘制背景场
    pc = ax.pcolormesh(glorys_lon_filtered, glorys_lat_filtered, glorys_variable_filtered, cmap='coolwarm', shading='auto', alpha=0.5)
    cbar = plt.colorbar(pc, ax=ax, orientation='horizontal', fraction=0.046, pad=0.04)
    if variable == 'vorticity':
        cbar.set_label(r'$\zeta/f$', fontsize=20)
        pc.set_clim(-0.8,0.8)
    elif variable == 'thetao':
        cbar.set_label('Temperature (°C)', fontsize=20)
    elif variable == 'so':
        cbar.set_label('Salinity (psu)', fontsize=20)
    elif variable == 'u':
        cbar.set_label('Zonal Velocity (m/s)', fontsize=20)
    elif variable == 'v':
        cbar.set_label('Meridional Velocity (m/s)', fontsize=20)
    elif variable == 'ssh':
        cbar.set_label('Sea Surface Height (m)', fontsize=20)
    else:
        cbar.set_label(variable, fontsize=20)
    cbar.ax.tick_params(labelsize=14)

    # 绘制Argo浮标数据
    if not needed_data.empty:
        if deep_argo:
            sc = ax.scatter(needed_data['Longitude'], needed_data['Latitude'], c=needed_data['DO_mol_kg'], cmap = 'bwr', s=180,
                            vmin=150, vmax=240, edgecolors='black', linewidths=0.5, label='Argo with max DO under 500m', zorder=5)
            cbar2 = plt.colorbar(sc, ax=ax, orientation='horizontal', fraction=0.046, pad=0.04)
            cbar2.set_label('DO/μmol·kg⁻¹', fontsize=20)
            cbar2.ax.tick_params(labelsize=14)
        else:
            ax.scatter(needed_data['Longitude'], needed_data['Latitude'], color='blue', s=180, label='Argo', zorder=5)
        for idx, row in needed_data.iterrows():
            ax.text(row['Longitude'], row['Latitude'], f"{int(row['Depth_m'])}", fontsize=7, fontweight='bold', ha='center', va='center', color='black', zorder=6)
    else:
        print(f"No Argo data available for eddy {ds_names}{no} at the specified index {needed_idx}.")

    # 绘制当前时刻涡旋
    circle = plt.Circle((center_lon[needed_idx], center_lat[needed_idx]),
                                radius[needed_idx] / 111320.0,
                                color='r', fill=False, linestyle='--', alpha=0.2, linewidth=1, label='Effective Radius')
    ax.add_patch(circle)

    ax.scatter(center_lon[needed_idx], center_lat[needed_idx], color='black', s=20, label='Eddy Center', zorder=5)
    ax.plot(contour_lon[needed_idx], contour_lat[needed_idx], color=colors, linewidth=1, alpha=0.5, label='Effective Contour')

    # 绘制 y = kx + b 直线
    if k is not None and b is not None:
        # 获取当前x轴的范围
        x_min, x_max = ax.get_xlim()
        # 根据 y = kx + b 计算y的对应范围
        line_x = np.array([x_min, x_max])
        line_y = k * line_x + b
        ax.plot(line_x, line_y, color='purple', linestyle='-', linewidth=2, label=f'Profile Line: y={k:.2f}x+{b:.2f}')

    ax.legend(fontsize=18)
    ax.set_xlim(glorys_lon_min, glorys_lon_max)
    ax.set_ylim(glorys_lat_min, glorys_lat_max)
    ax.set_aspect('equal')

    # 保存和显示图片
    if save_fig:
        output_dir = "plot_track_area_surface_glorys"
        os.makedirs(output_dir, exist_ok=True)
        
        base_filename = f"{ds_names}{no}_surface_{variable}_{dates[needed_idx].strftime('%Y%m%d')}.png"
        plt.savefig(os.path.join(output_dir, base_filename), dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {os.path.join(output_dir, base_filename)}")

    if show_fig:
        plt.show()

    plt.close(fig)  # 关闭图像，释放内存

def get_track_area_glorys(DS: list, no: int, needed_idx: int, variables: list = ['thetao'], surface: bool = False):
    '''
    获取指定涡旋在特定时间点周围的 GLORYS 数据。

    该函数会根据涡旋轮廓确定一个矩形区域，并从相应的 GLORYS 文件中
    提取此区域内的一个或多个物理变量。

    参数:
        DS (list): 包含涡旋轨迹信息的数据集。
        no (int): 涡旋的唯一编号。
        needed_idx (int): 涡旋轨迹的时间点索引。
        surface (bool): 是否只提取表层数据，默认为 False，提取2000m以内的数据。
        variables (list): 需要提取的变量列表，默认为 ['thetao']，可选'salinity', 'u', 'v', 'ssh'。

    返回:
        一个元组，包含筛选后的经度、纬度、深度数组，以及一个存储了所有请求变量数据的字典。
    '''
    wanted_track = find_track(DS, no)
    num, time, center_lon, center_lat, max_lon, max_lat, contour_lon, contour_lat, radius, speed_contour_lon, speed_contour_lat = zip(*wanted_track)

    glorys_filepaths_dict = find_track_glorys_filepath(DS, no)

    contour_lon_filtered = np.ma.masked_equal(contour_lon, 180.0)
    contour_lat_filtered = np.ma.masked_equal(contour_lat, 0.0)

    glorys_lon_min = np.min(contour_lon_filtered) - 0.5
    glorys_lon_max = np.max(contour_lon_filtered) + 0.5
    glorys_lat_min = np.min(contour_lat_filtered) - 0.5
    glorys_lat_max = np.max(contour_lat_filtered) + 0.5

    needed_glorys_data = Dataset(list(glorys_filepaths_dict.values())[needed_idx], 'r')
    glorys_lon = needed_glorys_data.variables['longitude'][:]
    glorys_lat = needed_glorys_data.variables['latitude'][:]
    glorys_depth = needed_glorys_data.variables['depth'][:]

    glorys_lon_mask = (glorys_lon >= glorys_lon_min) & (glorys_lon <= glorys_lon_max)
    glorys_lat_mask = (glorys_lat >= glorys_lat_min) & (glorys_lat <= glorys_lat_max)
    if surface:
        glorys_depth_mask = np.zeros_like(glorys_depth, dtype=bool)
        if glorys_depth.size > 0:
            glorys_depth_mask[0] = True
    else:
        glorys_depth_mask = (glorys_depth >= 0) & (glorys_depth <= 2000)
        
    glorys_lon_filtered = glorys_lon[glorys_lon_mask]
    glorys_lat_filtered = glorys_lat[glorys_lat_mask]
    glorys_depth_filtered = glorys_depth[glorys_depth_mask]

    # 存储多个变量的字典
    glorys_variables_filtered = {}
    for var in variables:
        if var == 'thetao':
            glorys_variables_filtered['thetao'] = needed_glorys_data.variables['thetao'][0, glorys_depth_mask, glorys_lat_mask, glorys_lon_mask]
            if surface:
                glorys_variables_filtered['thetao'] = glorys_variables_filtered['thetao'][0, :, :]
        elif var == 'salinity' or var == 'so':
            glorys_variables_filtered['salinity'] = needed_glorys_data.variables['so'][0, glorys_depth_mask, glorys_lat_mask, glorys_lon_mask]
            if surface:
                glorys_variables_filtered['salinity'] = glorys_variables_filtered['salinity'][0, :, :]
        elif var == 'u' or var == 'uo':
            glorys_variables_filtered['u'] = needed_glorys_data.variables['uo'][0, glorys_depth_mask, glorys_lat_mask, glorys_lon_mask]
            if surface:
                glorys_variables_filtered['u'] = glorys_variables_filtered['u'][0, :, :]
        elif var == 'v' or var == 'vo':
            glorys_variables_filtered['v'] = needed_glorys_data.variables['vo'][0, glorys_depth_mask, glorys_lat_mask, glorys_lon_mask]
            if surface:
                glorys_variables_filtered['v'] = glorys_variables_filtered['v'][0, :, :]
        elif var == 'ssh' or var == 'zos':
            glorys_variables_filtered['ssh'] = needed_glorys_data.variables['zos'][0, glorys_lat_mask, glorys_lon_mask]

    needed_glorys_data.close()
    
    return glorys_lon_filtered, glorys_lat_filtered, glorys_depth_filtered, glorys_variables_filtered

def calculate_vorticity(lon, lat, u, v):
    '''
    计算给定速度场的相对涡度 (zeta) 和科里奥利参数 (f)。
    该函数可智能处理2D/3D以及Masked Array输入。

    输入的速度场 u, v 的维度应为 (latitude, longitude) 或
    (depth, latitude, longitude)。如果输入为Masked Array，则输出也会是
    Masked Array，其中在mask边缘计算不准确的点会被自动mask掉。

    参数:
        lon (np.ndarray): 一维经度数组 (单位: 度)。
        lat (np.ndarray): 一维纬度数组 (单位: 度)。
        u (np.ndarray | np.ma.MaskedArray): Zonal 速度数组。
        v (np.ndarray | np.ma.MaskedArray): Meridional 速度数组。

    返回:
        tuple: 包含两个元素的元组 (zeta, f)。
               zeta (np.ndarray | np.ma.MaskedArray): 计算得到的相对涡度数组。
               f (np.ndarray | np.ma.MaskedArray): 计算得到的科里奥利参数数组。
    '''
    # --- 1. 输入校验和维度处理 ---
    if u.shape[-2:] != (len(lat), len(lon)) or u.shape != v.shape:
        raise ValueError("速度数组的最后两个维度必须与经纬度数组的长度匹配，且u,v数组形状需一致。")
    if u.ndim not in [2, 3]:
        raise ValueError("输入速度场必须是 2D 或 3D 数组。")

    # --- 2. 智能处理 Masked Array ---
    # 检查输入是否为 masked array
    is_masked_input = np.ma.is_masked(u)

    if is_masked_input:
        # 如果是，用 NaN 填充被 mask 的位置。梯度计算会正确传播NaN。
        u_data = u.filled(np.nan)
        v_data = v.filled(np.nan)
    else:
        # 如果不是，直接使用原始数据
        u_data = u
        v_data = v

    # --- 3. 统一升维处理 ---
    original_ndim = u.ndim
    if original_ndim == 2:
        u_proc = u_data[np.newaxis, :, :]
        v_proc = v_data[np.newaxis, :, :]
    else:
        u_proc = u_data
        v_proc = v_data

    # --- 4. 计算物理坐标间距 (dx, dy) 和科里奥利参数 (f) ---
    R_earth = 6371e3  # 地球半径 (米)
    Omega = 7.2921e-5  # 地球自转角速度 (弧度/秒)

    lon_rad, lat_rad = np.meshgrid(np.deg2rad(lon), np.deg2rad(lat))

    dy = R_earth * np.gradient(lat_rad, axis=0)
    dx = R_earth * np.cos(lat_rad) * np.gradient(lon_rad, axis=1)
    
    # 计算科里奥利参数 f
    # f 的形状应该与 (latitude, longitude) 匹配，然后根据需要广播到 (depth, latitude, longitude)
    f_2d = 2 * Omega * np.sin(lat_rad)

    # 广播 f_2d 到 u_proc 和 v_proc 的深度维度，以便最终输出的 f 形状与 zeta 保持一致
    if original_ndim == 3:
        f_proc = np.tile(f_2d[np.newaxis, :, :], (u_proc.shape[0], 1, 1))
    else: # original_ndim == 2
        f_proc = f_2d # f_2d 已经是 (latitude, longitude) 形状

    # --- 5. 核心计算逻辑 ---
    vorticity_list = []
    for k_slice_idx in range(u_proc.shape[0]):
        u_slice = u_proc[k_slice_idx, :, :]
        v_slice = v_proc[k_slice_idx, :, :]
        dvdx = np.gradient(v_slice, axis=1) / dx
        dudy = np.gradient(u_slice, axis=0) / dy
        vorticity_list.append(dvdx - dudy)
    vorticity_result = np.stack(vorticity_list, axis=0)
    
    # --- 6. 根据输入类型，决定最终输出 ---
    if is_masked_input:
        # 如果输入是 masked，将结果中的 NaN 转回为 mask
        zeta_final = np.ma.masked_invalid(vorticity_result, copy=False)
        f_final = np.ma.masked_invalid(f_proc, copy=False) # f 也可能是 masked array
    else:
        zeta_final = vorticity_result
        f_final = f_proc

    # 如果原始输入是 2D，则降维回去
    if original_ndim == 2:
        return zeta_final.squeeze(axis=0), f_final
    else:
        return zeta_final, f_final
    
def get_idx(DS: list, no: int, start_date: str, end_date: str) -> list:
    '''
    获取指定涡旋编号在给定时间范围内的索引列表。

    参数:
        DS (list): 涡旋轨迹数据集。
        no (int): 涡旋编号。
        start_date (str): 起始日期，格式为 'YYYY-MM-DD'。
        end_date (str): 结束日期，格式为 'YYYY-MM-DD'。

    返回:
        list: 涡旋编号在指定时间范围内的索引列表。
    '''
    wanted_track = find_track(DS, no)
    if not wanted_track:
        print(f"未找到涡旋 {no} 的轨迹数据。")
        return []

    start_date = pd.to_datetime(start_date)
    end_date = pd.to_datetime(end_date)

    idx_list = []
    for idx, track_point in enumerate(wanted_track):
        track_date = convert_date(track_point[1])
        if start_date <= track_date <= end_date:
            idx_list.append(idx)

    return idx_list

def plot_vertical_glorys(DS: list, no: int, needed_idx: int, k: float, b: float, variable: str = 'vorticity',
                         show_fig: bool = False, save_fig: bool = False, xmin: float = None, xmax: float = None,
                         ymin: float = None, ymax: float = None):
    '''
    绘制指定涡旋在特定时刻的垂直剖面图 (y = kx + b)。

    该函数会提取 GLORYS 数据的三维场，并在由 y = kx + b 定义的直线上进行插值，
    生成该变量的二维垂直剖面图像。

    参数:
        DS (list): 包含涡旋轨迹信息的数据集。
        no (int): 涡旋的唯一编号。
        needed_idx (int): 涡旋轨迹的时间点索引。
        k (float): 直线方程 y = kx + b 中的斜率。
        b (float): 直线方程 y = kx + b 中的截距。
        variable (str): 作为背景场绘制的GLORYS物理变量。默认为 'vorticity'。
        show_fig (bool): 是否在运行时显示生成的图像。默认为 False。
        save_fig (bool): 是否将生成的图像保存为 False。
        xmin (float): x轴最小值，默认为 None，表示自动计算。
        xmax (float): x轴最大值，默认为 None，表示自动计算。
        ymin (float): y轴最小值，默认为 None，表示自动计算。
        ymax (float): y轴最大值，默认为 None，表示自动计算。
    '''
    if k is None or b is None:
        raise ValueError("k 和 b 必须提供以绘制垂直剖面。")

    R_earth = 6371e3 # 地球半径 (米)

    wanted_track = find_track(DS, no)
    num, time, center_lon, center_lat, max_lon, max_lat, contour_lon, contour_lat, radius, speed_contour_lon, speed_contour_lat = zip(*wanted_track)
    dates = convert_date(time) if time else None

    # 获取涡旋中心和有效半径（度为单位）
    current_center_lon = center_lon[needed_idx]
    current_center_lat = center_lat[needed_idx]
    effective_radius_deg = radius[needed_idx] / 111320.0 # 将米转换为度

    # 获取区域边界
    contour_lon_filtered = np.ma.masked_equal(contour_lon, 180.0)
    contour_lat_filtered = np.ma.masked_equal(contour_lat, 0.0)

    glorys_lon_min = np.min(contour_lon_filtered) - 0.5
    glorys_lon_max = np.max(contour_lon_filtered) + 0.5
    glorys_lat_min = np.min(contour_lat_filtered) - 0.5
    glorys_lat_max = np.max(contour_lat_filtered) + 0.5
    
    # 获取三维背景场数据
    if variable == 'vorticity':
        glorys_lon_raw, glorys_lat_raw, glorys_depth_raw, glorys_variables_raw = \
            get_track_area_glorys(DS, no, needed_idx, variables=['u', 'v'], surface=False)
        
        u_data_for_vorticity = glorys_variables_raw['u']
        v_data_for_vorticity = glorys_variables_raw['v']
        if u_data_for_vorticity.ndim == 2:
            u_data_for_vorticity = u_data_for_vorticity[np.newaxis, :, :]
            v_data_for_vorticity = v_data_for_vorticity[np.newaxis, :, :]

        zeta_3d, f_3d = calculate_vorticity(glorys_lon_raw, glorys_lat_raw,
                                            u_data_for_vorticity, v_data_for_vorticity)
        glorys_variable_3d = zeta_3d / f_3d
        cbar_label = r'$\zeta/f$'
        cmap = 'coolwarm'
        clim = (-0.3, 0.3)
    else:
        glorys_lon_raw, glorys_lat_raw, glorys_depth_raw, glorys_variables_raw = \
            get_track_area_glorys(DS, no, needed_idx, variables=[variable], surface=False)
        
        if variable == 'salinity':
            glorys_variable_3d = glorys_variables_raw['salinity']
        elif variable == 'u':
            glorys_variable_3d = glorys_variables_raw['u']
        elif variable == 'v':
            glorys_variable_3d = glorys_variables_raw['v']
        elif variable == 'thetao':
            glorys_variable_3d = glorys_variables_raw['thetao']
        elif variable == 'ssh' or variable == 'zos':
            raise ValueError(f"Variable '{variable}' is a 2D surface variable. Vertical profile is not applicable.")
        else:
            glorys_variable_3d = glorys_variables_raw.get(variable)
            if glorys_variable_3d is None:
                 raise ValueError(f"Unknown or unsupported variable for vertical profile: '{variable}'.")

        if glorys_variable_3d.ndim == 2 and len(glorys_depth_raw) == 1:
            glorys_variable_3d = glorys_variable_3d[np.newaxis, :, :]
        elif glorys_variable_3d.ndim != 3:
             raise ValueError(f"Expected 3D data for variable '{variable}' but got {glorys_variable_3d.ndim}D data.")

        valid_values = glorys_variable_3d[~np.isnan(glorys_variable_3d)]
        if valid_values.size > 0:
            clim = (np.min(valid_values), np.max(valid_values))
        else:
            clim = (0, 1)

        if variable == 'thetao':
            cbar_label = 'Temperature (°C)'
            cmap = 'rainbow'
        elif variable == 'salinity':
            cbar_label = 'Salinity (psu)'
            cmap = 'viridis'
        elif variable == 'u' or variable == 'v':
            cbar_label = 'Velocity (m/s)'
            cmap = 'RdBu_r'
            if valid_values.size > 0:
                max_abs_val = np.max(np.abs(valid_values))
                clim = (-max_abs_val, max_abs_val)
            else:
                clim = (-1, 1)
        else:
            cbar_label = variable
            cmap = 'viridis'

    if glorys_variable_3d.size == 0 or np.all(np.isnan(glorys_variable_3d)):
        print(f"Warning: GLORYS data for variable '{variable}' is empty or all NaN for Track {no} at index {needed_idx}. Cannot plot vertical profile.")
        plt.close()
        return

    # --- 定义剖面线上的点 ---
    num_points = 500 # 增加点数，使交点查找更精确，并保证连续性
    
    if k == 0: # 水平线 (y = b)
        profile_lons = np.linspace(glorys_lon_min, glorys_lon_max, num_points)
        profile_lats = np.full_like(profile_lons, b)
    elif np.isinf(k) or np.abs(k) > 1e10: # 近似垂直线 (x = const)
        profile_lats = np.linspace(glorys_lat_min, glorys_lat_max, num_points)
        profile_lons = np.full_like(profile_lats, (glorys_lon_min + glorys_lon_max)/2.0)
    else: # 斜线
        profile_lons_temp = np.linspace(glorys_lon_min, glorys_lon_max, num_points)
        profile_lats_temp = k * profile_lons_temp + b
        
        valid_mask = (profile_lats_temp >= glorys_lat_min) & (profile_lats_temp <= glorys_lat_max)
        profile_lons = profile_lons_temp[valid_mask]
        profile_lats = profile_lats_temp[valid_mask]

        if len(profile_lons) < 2: # 如果经度范围不足以生成足够点，尝试以纬度为主
            profile_lats_temp = np.linspace(glorys_lat_min, glorys_lat_max, num_points)
            profile_lons_temp = (profile_lats_temp - b) / k
            valid_mask = (profile_lons_temp >= glorys_lon_min) & (profile_lons_temp <= glorys_lon_max)
            profile_lons = profile_lons_temp[valid_mask]
            profile_lats = profile_lats_temp[valid_mask]

    if len(profile_lons) < 2:
        print(f"Warning: The specified line (y={k}x+{b}) does not sufficiently intersect the GLORYS data area. Cannot plot vertical profile.")
        plt.close()
        return

    # === 修正 X 轴的物理距离计算 ===
    # 将经纬度转换为弧度
    profile_lons_rad = np.deg2rad(profile_lons)
    profile_lats_rad = np.deg2rad(profile_lats)

    # 计算每小段的经度差和纬度差 (弧度)
    dlon_rad_segments = np.diff(profile_lons_rad)
    dlat_rad_segments = np.diff(profile_lats_rad)

    # 计算每小段的平均纬度 (用于计算cos因子)
    mid_lats_rad = np.deg2rad((profile_lats[:-1] + profile_lats[1:]) / 2)

    # 将弧度差转换为米 (近似笛卡尔坐标)
    dx_meters = R_earth * np.cos(mid_lats_rad) * dlon_rad_segments
    dy_meters = R_earth * dlat_rad_segments

    # 计算每小段的斜边距离 (米)
    distance_segments_meters = np.sqrt(dx_meters**2 + dy_meters**2)

    # 累积距离 (公里)
    x_coords_raw = np.insert(np.cumsum(distance_segments_meters), 0, 0) / 1000.0
    # ==========================================================

    # === 1. 计算涡旋中心在剖面上的投影点 ===
    if k == 0: # 水平线 y=b
        proj_center_lon = current_center_lon
        proj_center_lat = b
    elif np.isinf(k) or np.abs(k) > 1e10: # 近似垂直线 x=const
        proj_center_lon = (glorys_lon_min + glorys_lon_max)/2.0 # 或特定x值
        proj_center_lat = current_center_lat
    else: # 斜线 y = kx + b
        xp = (current_center_lon + k * current_center_lat - k * b) / (1 + k*k)
        yp = k * xp + b
        proj_center_lon = xp
        proj_center_lat = yp

    # 找到 (proj_center_lon, proj_center_lat) 在 profile_lons, profile_lats 序列中的索引
    # 先将 profile_lons, profile_lats 转为点对
    profile_points_xy = np.column_stack((profile_lons, profile_lats))
    proj_point_xy = np.array([proj_center_lon, proj_center_lat])

    # 计算投影点到剖面线上每个点的欧氏距离（在经纬度空间）
    distances_to_profile_points_lonlat = np.sqrt(
        (profile_points_xy[:, 0] - proj_point_xy[0])**2 +
        (profile_points_xy[:, 1] - proj_point_xy[1])**2
    )
    # 找到距离最近的点的索引
    closest_point_idx_on_profile = np.argmin(distances_to_profile_points_lonlat)
    
    # 涡旋中心投影点在 x_coords_raw 中的距离
    center_proj_dist_from_start = x_coords_raw[closest_point_idx_on_profile]

    # 重新中心化 X 轴
    x_coords_recenter = x_coords_raw - center_proj_dist_from_start
    # ==========================================================

    # --- 对三维数据进行插值以获取剖面数据 ---
    if isinstance(glorys_variable_3d, np.ma.MaskedArray):
        interp_data = glorys_variable_3d.filled(np.nan)
    else:
        interp_data = glorys_variable_3d

    points_grid = (glorys_depth_raw, glorys_lat_raw, glorys_lon_raw)
    interp_func = RegularGridInterpolator(points_grid, interp_data, method='linear', bounds_error=False, fill_value=np.nan)

    query_depths, query_lats = np.meshgrid(glorys_depth_raw, profile_lats, indexing='ij')
    _, query_lons = np.meshgrid(glorys_depth_raw, profile_lons, indexing='ij')

    xi_interp_points = np.vstack([query_depths.ravel(), query_lats.ravel(), query_lons.ravel()]).T

    interpolated_values_flat = interp_func(xi_interp_points)
    profile_variable_2d = interpolated_values_flat.reshape(len(glorys_depth_raw), len(profile_lons))

    if isinstance(glorys_variable_3d, np.ma.MaskedArray):
        profile_variable_2d = np.ma.masked_invalid(profile_variable_2d, copy=False)

    # --- 绘图 ---
    fig, ax = plt.subplots(figsize=(20, 15))
    
    callers_local_vars = inspect.currentframe().f_back.f_locals.items()
    ds_names = [var_name for var_name, var_val in callers_local_vars if var_val is DS]
    if ds_names:
        ds_names = ds_names[0].upper()
    else:
        ds_names = "UNKNOWN"

    prop_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    if ds_names == 'ACS' or ds_names == 'ACL':
        eddy_color = prop_colors[1]
    else:
        eddy_color = prop_colors[0]

    ax.set_title(f'Vertical Profile of {cbar_label} for Track {ds_names}{num[0]} on {dates[needed_idx].strftime('%Y-%m-%d')}', fontsize=20)
    
    ax.set_xlabel('Distance from Eddy Center Projection (km)', fontsize=18)
    ax.set_ylabel('Depth (m)', fontsize=18)
    ax.tick_params(labelsize=14)

    X_mesh_for_plot, Y_mesh_for_plot = np.meshgrid(x_coords_recenter, glorys_depth_raw)
    
    pc = ax.pcolormesh(X_mesh_for_plot, Y_mesh_for_plot, profile_variable_2d, cmap=cmap, shading='auto', vmin=clim[0], vmax=clim[1])
    cbar = plt.colorbar(pc, ax=ax, orientation='vertical', fraction=0.046, pad=0.04)
    cbar.set_label(cbar_label, fontsize=18)
    cbar.ax.tick_params(labelsize=14)

    ax.set_ylim(glorys_depth_raw.max(), glorys_depth_raw.min())

    # === 绘制涡旋中心投影的垂直线 ===
    ax.axvline(0, color='black', linestyle='--', linewidth=2, label='Eddy Center Projection')

    # === 绘制有效半径 (circle) 的垂直线 ===
    A_quad = 1 + k*k # 避免与A参数混淆
    B_quad = -2*current_center_lon + 2*k*(b - current_center_lat)
    C_quad = current_center_lon**2 + (b - current_center_lat)**2 - effective_radius_deg**2
    
    discriminant = B_quad**2 - 4*A_quad*C_quad
    
    circle_intersections_xy = []
    if discriminant >= 0:
        x1 = (-B_quad + np.sqrt(discriminant)) / (2*A_quad)
        y1 = k * x1 + b
        circle_intersections_xy.append((x1, y1))
        if discriminant > 0:
            x2 = (-B_quad - np.sqrt(discriminant)) / (2*A_quad)
            y2 = k * x2 + b
            circle_intersections_xy.append((x2, y2))
            
    circle_proj_distances = []
    for int_lon, int_lat in circle_intersections_xy:
        # 找到交点在 profile_lons, profile_lats 序列中的索引
        int_point_xy = np.array([int_lon, int_lat])
        distances_to_profile_points_from_int = np.sqrt(
            (profile_points_xy[:, 0] - int_point_xy[0])**2 + 
            (profile_points_xy[:, 1] - int_point_xy[1])**2
        )
        closest_profile_idx = np.argmin(distances_to_profile_points_from_int)
        
        dist_from_start = x_coords_raw[closest_profile_idx]
        relative_dist = dist_from_start - center_proj_dist_from_start
        
        # 仅当交点接近剖面线时才绘制
        if distances_to_profile_points_from_int[closest_profile_idx] < effective_radius_deg * 0.1: # 小于半径的10%作为阈值
             circle_proj_distances.append(relative_dist)
            
    for dist in sorted(list(set(circle_proj_distances))): # 去重并排序
        ax.axvline(dist, color='r', linestyle='--', linewidth=1, label='Effective Radius Projection' if dist == sorted(list(set(circle_proj_distances)))[0] else "")


    # === 绘制有效轮廓 (contour) 的垂直线 ===
    current_contour_lon = np.array(contour_lon)[needed_idx]
    current_contour_lat = np.array(contour_lat)[needed_idx]

    valid_contour_mask = (current_contour_lon != 180.0) & (current_contour_lat != 0.0)
    current_contour_lon_valid = current_contour_lon[valid_contour_mask]
    current_contour_lat_valid = current_contour_lat[valid_contour_mask]

    if len(current_contour_lon_valid) > 2:
        contour_intersections_xy = find_polygon_line_intersections(
            current_contour_lon_valid, current_contour_lat_valid,
            profile_lons, profile_lats
        )
        
        contour_proj_distances = []
        for int_lon, int_lat in contour_intersections_xy:
            int_point_xy = np.array([int_lon, int_lat])
            distances_to_profile_points_from_int = np.sqrt(
                (profile_points_xy[:, 0] - int_point_xy[0])**2 + 
                (profile_points_xy[:, 1] - int_point_xy[1])**2
            )
            closest_profile_idx = np.argmin(distances_to_profile_points_from_int)
            
            dist_from_start = x_coords_raw[closest_profile_idx]
            relative_dist = dist_from_start - center_proj_dist_from_start
            
            # 检查交点是否在剖面线的有效X范围（防止绘制图外交点）
            if np.min(x_coords_recenter) <= relative_dist <= np.max(x_coords_recenter):
                contour_proj_distances.append(relative_dist)
        
        for dist in sorted(list(set(contour_proj_distances))): # 去重并排序
            ax.axvline(dist, color=eddy_color, linestyle=':', linewidth=1, label='Effective Contour Projection' if dist == sorted(list(set(contour_proj_distances)))[0] else "")
    else:
        print(f"Warning: Not enough valid points to define effective contour for Track {no} at index {needed_idx}.")

    # 统一图例，避免重复标签
    handles, labels = ax.get_legend_handles_labels()
    unique_labels = {}
    for handle, label in zip(handles, labels):
        unique_labels[label] = handle
    ax.legend(unique_labels.values(), unique_labels.keys(), fontsize=18)

    if xmin is not None and xmax is not None:
        ax.set_xlim(xmin, xmax)
    if ymin is not None and ymax is not None:
        ax.set_ylim(ymax, ymin)

    # 保存和显示图片
    if save_fig:
        output_dir = "plot_vertical_glorys"
        os.makedirs(output_dir, exist_ok=True)
        base_filename = f"{ds_names}{no}_vertical_{variable}_{dates[needed_idx].strftime('%Y%m%d')}_k{k:.2f}b{b:.2f}.png"
        plt.savefig(os.path.join(output_dir, base_filename), dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {os.path.join(output_dir, base_filename)}")

    if show_fig:
        plt.show()

    plt.close(fig)

# 帮助函数：判断三个点 (p, q, r) 的方向（共线，顺时针，逆时针）
def _orientation(p, q, r):
    val = (q[1] - p[1]) * (r[0] - q[0]) - (q[0] - p[0]) * (r[1] - q[1])
    if val == 0: return 0  # Collinear
    return 1 if val > 0 else 2 # Clockwise or Counterclockwise (1 for clockwise, 2 for counterclockwise)

# 帮助函数：判断点 q 是否在线段 pr 上
def _on_segment(p, q, r):
    return (q[0] <= max(p[0], r[0]) and q[0] >= min(p[0], r[0]) and
            q[1] <= max(p[1], r[1]) and q[1] >= min(p[1], r[1]))

# 核心帮助函数：计算两条线段的交点
def _line_segment_intersect(p1, q1, p2, q2):
    """
    Finds the intersection point of two line segments (p1, q1) and (p2, q2).
    Returns the intersection point as (x, y) if they intersect, otherwise None.
    Handles collinear cases where segments overlap.
    """
    o1 = _orientation(p1, q1, p2)
    o2 = _orientation(p1, q1, q2)
    o3 = _orientation(p2, q2, p1)
    o4 = _orientation(p2, q2, q1)

    # General case: Non-collinear and non-parallel intersection
    if o1 != o2 and o3 != o4:
        denom = (p1[0] - q1[0]) * (p2[1] - q2[1]) - (p1[1] - q1[1]) * (p2[0] - q2[0])
        if denom == 0:
            return None 

        t = ((p1[0] - p2[0]) * (p2[1] - q2[1]) - (p1[1] - p2[1]) * (p2[0] - q2[0])) / denom
        u = -((p1[0] - q1[0]) * (p1[1] - p2[1]) - (p1[1] - q1[1]) * (p1[0] - p2[0])) / denom

        if (0 <= t <= 1) and (0 <= u <= 1):
            intersection_x = p1[0] + t * (q1[0] - p1[0])
            intersection_y = p1[1] + t * (q1[1] - p1[1])
            return (intersection_x, intersection_y)
        return None

    # Special Cases for collinear segments
    if o1 == 0 and _on_segment(p1, p2, q1): return p2
    if o2 == 0 and _on_segment(p1, q2, q1): return q2
    if o3 == 0 and _on_segment(p2, p1, q2): return p1
    if o4 == 0 and _on_segment(p2, q1, q2): return q2

    return None


# 寻找多边形和线段的交点
def find_polygon_line_intersections(polygon_lon, polygon_lat, line_lons, line_lats, tolerance=1e-6):
    """
    Finds intersection points of a line (defined by line_lons, line_lats)
    with a closed polygon (defined by polygon_lon, polygon_lat).

    Args:
        polygon_lon (array): Longitudes of polygon vertices.
        polygon_lat (array): Latitudes of polygon vertices.
        line_lons (array): Longitudes of line points.
        line_lats (array): Latitudes of line points.
        tolerance (float): Tolerance for checking duplicate intersection points.

    Returns:
        list: A list of (lon, lat) tuples for intersection points.
    """
    intersections = []
    
    polygon_points_closed = list(zip(polygon_lon, polygon_lat))
    if not (np.isclose(polygon_points_closed[0][0], polygon_points_closed[-1][0], atol=tolerance) and
            np.isclose(polygon_points_closed[0][1], polygon_points_closed[-1][1], atol=tolerance)):
        polygon_points_closed.append(polygon_points_closed[0])

    line_points = list(zip(line_lons, line_lats))

    for i in range(len(polygon_points_closed) - 1):
        poly_seg_p1 = polygon_points_closed[i]
        poly_seg_p2 = polygon_points_closed[i+1]
        
        for j in range(len(line_points) - 1):
            line_seg_p1 = line_points[j]
            line_seg_p2 = line_points[j+1]
            
            intersect_pt = _line_segment_intersect(poly_seg_p1, poly_seg_p2, line_seg_p1, line_seg_p2)
            if intersect_pt:
                is_duplicate = False
                for existing_pt in intersections:
                    if np.isclose(existing_pt[0], intersect_pt[0], atol=tolerance) and \
                       np.isclose(existing_pt[1], intersect_pt[1], atol=tolerance):
                        is_duplicate = True
                        break
                if not is_duplicate:
                    intersections.append(intersect_pt)
    
    if len(intersections) > 0:
        if np.abs(line_lons[0] - line_lons[-1]) > np.abs(line_lats[0] - line_lats[-1]):
            intersections.sort(key=lambda p: p[0])
        else:
            intersections.sort(key=lambda p: p[1])
            
    return intersections