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

lonmin,lonmax=140-2.5, 180+2.5
latmin,latmax=28-2.5, 40+2.5

argo_data = pd.read_parquet("Argo2014.parquet")
argo_data = argo_data.drop(columns=['Salinity_psu', 'Oxygen_flag', 'Oxygen_flag2', 'Datasets_number', 'Cycle_number', 'Float_serial_no'])

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
    
    筛选方式：先寻找指定编号涡旋出现日期，由日期筛选浮标数据，再判断浮标是否在涡旋contour或是effective radius内
    '''
    circle_enlargement_factor = 1.2  # 筛选过程中涡旋半径放大倍数

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
    绘制涡旋轨迹与浮标数据
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

    # For each observation in needed_data, find the matching track point (by date)
    labeled = False
    for idx, row in needed_data.iterrows():
        obs_date = pd.Timestamp(year=int(row['Year']),
                                month=int(row['Month']),
                                day=int(row['Day']))
        matching = np.where(pd.to_datetime(dates) == obs_date)[0]
        if matching.size > 0:
            i = matching[0]
            # Draw the radius as a circle (convert meters to degrees approximately)
            if labeled == False:
                if ds_names == 'ACS' or ds_names == 'ACL':
                    circle = plt.Circle((center_lon[i], center_lat[i]),
                                radius[i] / 111320,
                                color='r', fill=False, linestyle='--', alpha=0.2, linewidth=1, label='Effective Radius')
                else:
                    circle = plt.Circle((center_lon[i], center_lat[i]),
                                radius[i] / 111320,
                                color='purple', fill=False, linestyle='--', alpha=0.2, linewidth=1, label='Effective Radius')        
            else:
                if ds_names == 'ACS' or ds_names == 'ACL':
                    circle = plt.Circle((center_lon[i], center_lat[i]),
                                    radius[i] / 111320,
                                    color='r', fill=False, linestyle='--', alpha=0.2, linewidth=1)
                else:
                    circle = plt.Circle((center_lon[i], center_lat[i]),
                                    radius[i] / 111320,
                                    color='purple', fill=False, linestyle='--', alpha=0.2, linewidth=1)
            ax.add_patch(circle)
            # Mark the center point
            # ax.plot(center_lon[i], center_lat[i], marker='o', color='black', markersize=10)
            # Plot the contour line
            if labeled == False:
                ax.plot(contour_lon[i], contour_lat[i], color=colors, linewidth=1, alpha=0.5, label='Effective Contour')
                labeled = True
            else:
                ax.plot(contour_lon[i], contour_lat[i], color=colors, linewidth=1, alpha=0.5)
            # Plot date
            dates_obs = pd.to_datetime(needed_data[['Year', 'Month', 'Day']].astype(int))
            if obs_date == dates_obs.min() or obs_date == dates_obs.max():
                ax.text(center_lon[i], center_lat[i], obs_date.strftime('%Y-%m-%d'), fontsize=20, color='black')

    # Plot each observation point from needed_data
    ax.scatter(needed_data['Longitude'], needed_data['Latitude'], color='black', s=50, label='Floats')
    
    x_margin = (max(center_lon) - min(center_lon)) * 0.1
    y_margin = (max(center_lat) - min(center_lat)) * 0.1
    ax.set_xlim(min(center_lon) - x_margin, max(center_lon) + x_margin)
    ax.set_ylim(min(center_lat) - y_margin, max(center_lat) + y_margin)

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