import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon
import geopandas as gpd
import inspect
from netCDF4 import Dataset
import os
import pickle

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

def plot_vertical(DS: list, no: int, show_fig: bool = False, save_fig: bool = False, color_mode: str = 'distance'):
    '''
    根据涡旋轨迹和浮标数据，绘制浮标的垂直剖面图。

    参数:
        DS (list): 涡旋轨迹数据集。
        no (int): 涡旋编号。
        show_fig (bool): 是否显示图片，默认False。
        save_fig (bool): 是否保存图片，默认False。
        color_mode (str): 颜色模式，可选'distance'（按距离着色）或'time'（按时间着色）。

    功能:
        对每个浮标平台，按Profile_number分组，绘制DO_mol_kg随深度变化的曲线。
        曲线颜色可根据浮标与涡旋中心的相对距离或采样时间变化。
        支持图片保存与显示。
    '''
    wanted_track = find_track(DS, no)
    argo_data_filtered = filtered_float_data(DS, no)

    callers_local_vars = inspect.currentframe().f_back.f_locals.items()     #加入会导致调试卡顿
    ds_names = [var_name for var_name, var_val in callers_local_vars if var_val is DS]
    if ds_names:
        ds_names = ds_names[0].upper()
    else:
        raise ValueError("UNKNOWN VARIABLE")
    
    for idx,platform in argo_data_filtered.groupby("Platform_number"):
        idxx_values = platform['Profile_number'].unique()
        min_idxx = min(idxx_values)
        max_idxx = max(idxx_values)
        

        fig = plt.figure(figsize=(30, 20))
        cmap = plt.cm.coolwarm
        for idxx,rows in platform.groupby("Profile_number"):
            if color_mode == 'distance':
                needed_data = rows.iloc[0]
                date = pd.Timestamp(year=int(rows.iloc[0]['Year']),
                                    month=int(rows.iloc[0]['Month']),
                                    day=int(rows.iloc[0]['Day']))
                track_date = convert_date([t[1] for t in wanted_track])
                idx_track = np.where(track_date == date)[0]
                
                if len(idx_track) > 0:
                    idx_track = idx_track[0]
                    center_lon, center_lat, radius = wanted_track[idx_track][2], wanted_track[idx_track][3], wanted_track[idx_track][8]
                    rel_x = (needed_data['Longitude'] - center_lon) / (radius / 111320)
                    rel_y = (needed_data['Latitude'] - center_lat) / (radius / 111320)
                    distance = np.sqrt(rel_x**2 + rel_y**2)
                    normalized_distance = (distance - 0) / (1 - 0 + 1e-8)
                    color = cmap(1-normalized_distance)
                else:
                    print(f"Date {date} not found in track data.")
            elif color_mode == 'time':
                normalized_idxx = (idxx - min_idxx) / (max_idxx - min_idxx)
                color = cmap(normalized_idxx)
            obs_date=pd.Timestamp(year=int(rows.iloc[0]['Year']), month=int(rows.iloc[0]['Month']), day=int(rows.iloc[0]['Day']))
            plt.plot(rows["DO_mol_kg"], rows["Depth_m"], label=obs_date.strftime("%Y-%m-%d"), color=color, alpha=0.7)
            
        date_start = pd.Timestamp(year=int(platform.iloc[0]['Year']),
                                    month=int(platform.iloc[0]['Month']),
                                    day=int(platform.iloc[0]['Day']))
        date_end = pd.Timestamp(year=int(platform.iloc[-1]['Year']),
                                month=int(platform.iloc[-1]['Month']),
                                day=int(platform.iloc[-1]['Day']))
        
        plt.title(f"{ds_names}{no}, Platform Number: {int(idx)}, {date_start.date()}~{date_end.date()}", fontsize=20)
        plt.xlabel("DO_mol/kg", fontsize=20)
        plt.ylabel("Depth/m", fontsize=20)
        plt.tick_params(axis='both', which='major', labelsize=16)
        plt.gca().invert_yaxis()
        # plt.legend()

        # 保存图片
        if save_fig:
            output_dir = "plot_vertical"
            os.makedirs(output_dir, exist_ok=True)
            plt.savefig(os.path.join(output_dir, f"vertical_profile_platform_{idx}.png"), dpi=300, bbox_inches='tight')

        if show_fig:
            plt.show()
        plt.close(fig)
        
def plot_relative_position(DS: list, no: int, show_fig: bool = False, save_fig: bool = False, color_mode: str = 'distance'):
    '''
    根据涡旋轨迹和浮标数据，绘制浮标在单位圆涡旋中的相对位置分布图。

    参数:
        DS (list): 涡旋轨迹数据集。
        no (int): 涡旋编号。
        show_fig (bool): 是否显示图片，默认False。
        save_fig (bool): 是否保存图片，默认False。
        color_mode (str): 颜色模式，可选'distance'（按距离着色）或'time'（按时间着色）。

    功能:
        对每个浮标平台，按Profile_number分组，计算每个观测点相对于涡旋中心和半径的归一化位置（单位圆内）。
        观测点颜色可根据与涡旋中心的相对距离或采样时间变化。
        绘制单位圆表示涡旋边界，中心点为涡旋中心。
        支持图片保存与显示。
    '''
    wanted_track = find_track(DS, no)
    argo_data_filtered = filtered_float_data(DS, no)

    callers_local_vars = inspect.currentframe().f_back.f_locals.items()     #加入会导致调试卡顿
    ds_names = [var_name for var_name, var_val in callers_local_vars if var_val is DS]
    if ds_names:
        ds_names = ds_names[0].upper()
    else:
        raise ValueError("UNKNOWN VARIABLE")

    for idx, platform in argo_data_filtered.groupby("Platform_number"):
                
        needed_data = platform.groupby("Profile_number").apply(lambda group: group.iloc[0])
        needed_data.index.name = None

        dates = pd.to_datetime({'year': needed_data["Year"],
                                'month': needed_data["Month"],
                                'day': needed_data["Day"]})
        # 计算每个dates在wanted_track中的idx_track
        track_dates = convert_date([t[1] for t in wanted_track])
        idx_track_list = []
        for d in dates:
            # 找到与d匹配的track日期索引
            matches = np.where(track_dates == d)[0]
            if len(matches) > 0:
                idx_track_list.append(matches[0])
            else:
                idx_track_list.append(None)  # 若无匹配则填None

        # 构建needed_track_data，每行为中心点经度，中心点纬度，半径
        needed_track_data = []
        for idx_track in idx_track_list:
            if idx_track is not None:
                needed_track_data.append([
                    wanted_track[idx_track][2],  # center_lon
                    wanted_track[idx_track][3],  # center_lat
                    wanted_track[idx_track][8]   # radius
                ])
            else:
                needed_track_data.append([None, None, None]) 

        # 绘制每个浮标相对单位圆涡旋的位置
        fig, ax = plt.subplots(figsize=(30, 20))
        date_start = pd.Timestamp(year=int(needed_data.iloc[0]['Year']),
                                  month=int(needed_data.iloc[0]['Month']),
                                  day=int(needed_data.iloc[0]['Day']))
        date_end = pd.Timestamp(year=int(needed_data.iloc[-1]['Year']),
                                month=int(needed_data.iloc[-1]['Month']),
                                day=int(needed_data.iloc[-1]['Day']))
        
        ax.set_title(f"{ds_names}{no}, Platform Number: {int(idx)}, {date_start.date()}~{date_end.date()}, Total: {len(needed_data)}", fontsize=20)
        ax.set_xlabel('Longitude', fontsize=20)
        ax.set_ylabel('Latitude', fontsize=20)
        cmap = plt.cm.coolwarm  # 使用渐变色

        for i, (index, row) in enumerate(needed_data.iterrows()):
            center_lon, center_lat, radius = needed_track_data[i]
            if None in (center_lon, center_lat, radius):
                continue  # 跳过没有匹配track的点
            # 计算相对位置（单位圆）
            rel_x = (row['Longitude'] - center_lon) / (radius / 111320)
            rel_y = (row['Latitude'] - center_lat) / (radius / 111320)
            # 颜色映射
            if color_mode == 'distance':
                distance = np.sqrt(rel_x**2 + rel_y**2)
                normalized_distance = (distance - 0) / (1 - 0 + 1e-8)
                color = cmap(1-normalized_distance)
            elif color_mode == 'time':
                normalized_idx = (index - needed_data.index.min()) / (needed_data.index.max() - needed_data.index.min() + 1e-8)
                color = cmap(normalized_idx)
            ax.scatter(rel_x, rel_y, color=color, s=300, label=pd.Timestamp(year=int(row['Year']), month=int(row['Month']), day=int(row['Day'])).strftime("%Y-%m-%d"))
            ax.plot(0, 0, marker='x', color='black', markersize=16, markeredgewidth=3, label='Eddy Center' if i == 0 else "")
            # 绘制顺序
            # ax.text(rel_x, rel_y, pd.Timestamp(year=int(row['Year']), month=int(row['Month']), day=int(row['Day'])).strftime("%Y-%m-%d"), fontsize=12, color='black', ha='center', va='center')
            ax.text(rel_x, rel_y, i+1, weight='bold', fontsize=7, color='black', ha='center', va='center')
            

        # 绘制单位圆
        circle = plt.Circle((0, 0), 1, color='gray', fill=False, linestyle='--', linewidth=2, label='Unit Eddy')
        ax.add_patch(circle)
        ax.set_aspect('equal')
        # ax.legend(fontsize=14)

        # 保存图片
        if save_fig:
            output_dir = "plot_relative_position"
            os.makedirs(output_dir, exist_ok=True)
            plt.savefig(os.path.join(output_dir, f"relative_position_platform_{idx}.png"), dpi=300, bbox_inches='tight')
        
        if show_fig:
            plt.show()
        plt.close(fig)