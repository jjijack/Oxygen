import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon
import geopandas as gpd
import inspect
from netCDF4 import Dataset
import os
from pathlib import Path
import pickle
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import glob
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import gaussian_filter
import copy
import gsw
from collections import defaultdict
import multiprocessing
import h5py
import time as tm
import shutil
from matplotlib.lines import Line2D
import dask.dataframe as dd
from dask.distributed import Client
from tqdm.auto import tqdm

lonmin,lonmax=140-2.5, 180+2.5
latmin,latmax=28-2.5, 40+2.5

argo_origin_path = Path("./Argo_origin")
tmp_parquet_path = Path("./Argo_data_tmp")
argo_path = Path("./Argo_data")

# # load 2014 data by default, can be changed with load_argo_data function
# default_argo_data_path = argo_path / 'Argo2014.parquet'
# try:
#     argo_data = pd.read_parquet(default_argo_data_path)
#     argo_data = argo_data.drop(columns=['Salinity_psu', 'Oxygen_flag', 'Oxygen_flag2', 'Datasets_number', 'Cycle_number', 'Float_serial_no'])
#     print("Old format Argo data loaded successfully.")
# except FileNotFoundError:
#     print(f"Default Argo data file not found at {default_argo_data_path}. Empty DataFrame created.")
#     argo_data = pd.DataFrame()
# except KeyError:
#     print("New format Argo data loaded successfully.")

circle_enlargement_factor = 1.2  # 筛选过程中涡旋半径放大倍数
Glorys_path = '../copernicus/GLORYS'

def load_meta_data(path, version = 3.2):
    '''
    加载meta数据，输出ACS, ACL, CS, CL四个数据集。
    默认版本为3.2
    '''
    if version == 3.2:
        ACS= Dataset(os.path.join(path, 'META3.2_DT_allsat_Anticyclonic_short_19930101_20220209.nc'))
        ACL= Dataset(os.path.join(path, 'META3.2_DT_allsat_Anticyclonic_long_19930101_20220209.nc'))
        CS= Dataset(os.path.join(path, 'META3.2_DT_allsat_Cyclonic_short_19930101_20220209.nc'))
        CL= Dataset(os.path.join(path, 'META3.2_DT_allsat_Cyclonic_long_19930101_20220209.nc'))
    elif version == 3.1:
        ACS=Dataset(os.path.join(path, 'META3.1exp_DT_allsat_Anticyclonic_short_19930101_20200307.nc'))
        ACL=Dataset(os.path.join(path, 'META3.1exp_DT_allsat_Anticyclonic_long_19930101_20200307.nc'))
        CS=Dataset(os.path.join(path, 'META3.1exp_DT_allsat_Cyclonic_short_19930101_20200307.nc'))
        CL=Dataset(os.path.join(path, 'META3.1exp_DT_allsat_Cyclonic_long_19930101_20200307.nc'))
    else:
        raise ValueError("Unsupported version. Please use 3.1 or 3.2.")

    return ACS, ACL, CS, CL


def area_limit(DS, latmin, latmax, lonmin, lonmax):
    '''
    限制区域范围
    
    输出：涡旋序号，时间，中心点经度，中心点纬度，最值点经度，最值点纬度，边界经度，边界纬度，半径，速度边界经度，速度边界纬度
    '''
    time = DS.variables['time'][:].data
    if np.issubdtype(time.dtype, np.floating):
        # 如果时间是浮点数，转换为整数
        time = np.round(time).astype(np.uint32)
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

def convert_mat_to_parquet(year: int, input_dir: str = './Argo_addFloat', output_dir: str = './Argo_data'):
    """
    读取指定年份的多个Argo .mat文件，将其合并，并保存为单个Parquet文件。

    参数:
        year (int): 需要处理的数据年份。
        input_dir (str, optional): 存放源 .mat 文件的目录。
                                   默认为 './Argo_addFloat'。
        output_dir (str, optional): 用于保存输出 .parquet 文件的目录。
                                    默认为 './Argo_data'。
    """
    # 将年份转换为字符串，用于构建路径
    year_str = str(year)
    print(f"Starting to process data for the year: {year_str}")

    # --- 1. 准备文件路径 ---
    # 构建源文件路径列表
    try:
        paths = [os.path.join(input_dir, f'Argo{year_str}_{month}.mat') for month in range(1, 13)]
        print("Generated paths for .mat files:")
        for p in paths:
            print(f"  - {p}")
    except Exception as e:
        print(f"Error creating file paths: {e}")
        return

    # --- 2. 读取并处理数据 ---
    # 定义列名
    columns = [
        "Year", "Month", "Day", "Longitude", "Latitude", "Depth_m", "DO_mol_kg", "Salinity_psu",
        "Temperature_degC", "Oxygen_flag", "Oxygen_flag2", "Profile_number", "Datasets_number",
        "Platform_number", "Cycle_number", "Float_serial_no"
    ]
    
    all_data = [] # 初始化空列表用于存储所有月份的 DataFrame

    print("\nReading .mat files...")
    for path in paths:
        if not os.path.exists(path):
            print(f"Warning: File not found, skipping: {path}")
            continue
        try:
            with h5py.File(path, 'r') as f:
                if 'do' in f:
                    do_data = f['do'][:].T  # 转置为行优先格式
                    df = pd.DataFrame(do_data, columns=columns)
                    all_data.append(df)
                    print(f"Successfully read and processed {path}")
                else:
                    print(f"Warning: 'do' dataset not found in {path}")
        except Exception as e:
            print(f"Error: Failed to read {path}. Reason: {e}")

    # --- 3. 合并并保存数据 ---
    if not all_data:
        print("\nNo data was loaded. Nothing to save.")
        return

    print("\nConcatenating all loaded data...")
    final_df = pd.concat(all_data, ignore_index=True)

    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)
    
    # 构建输出文件路径并保存
    output_path = os.path.join(output_dir, f"Argo{year_str}.parquet")
    try:
        final_df.to_parquet(output_path, index=False)
        print(f"\nSuccess! All data for {year_str} successfully saved to:")
        print(f"  -> {output_path}")
    except Exception as e:
        print(f"\nError: Failed to save Parquet file. Reason: {e}")

def parse_argo_txt_file(file_path: Path) -> pd.DataFrame | None:
    """
    解析单个扁平的、由制表符分隔的Argo表格文件。

    功能:
        此函数假设输入文件是一个标准的表格文件，其中第一行是列标题，
        后续行是数据，所有字段均由制表符 ('\t') 分隔。
        它会自动将指定的占位符（如 -999）转换成标准的 NaN。

    参数:
        file_path (Path): 单个Argo txt文件的路径。

    返回:
        pd.DataFrame | None: 如果文件成功解析且内容不为空，返回DataFrame；否则返回None。
    """
    try:
        # 定义需要被视作NaN的占位符列表
        missing_values = [-999, -999.0, -999.00, -999.000, -999.0000, -999.000000]
        
        # 在读取时使用 na_values 参数直接完成替换
        df = pd.read_csv(
            file_path,
            sep='\t',
            na_values=missing_values
        )
        
        # 如果文件只有表头或完全为空，df.empty会是True
        if df.empty:
            return None
        return df
    except Exception as e:
        # print(f"  [Warning] 读取文件 {file_path.name} 失败: {e}")
        return None

def worker_process_file(task_args: tuple[Path, Path]):
    """
    (Worker函数) 这是Dask的工作单元，负责处理单个文件。

    功能:
        接收一个包含输入和输出路径的元组，调用解析函数，
        并将结果保存为中间Parquet文件。

    参数:
        task_args (tuple[Path, Path]): 一个元组，包含:
            - 原始txt文件路径 (input_path)
            - 中间Parquet文件输出路径 (output_path)
    """
    input_path, output_path = task_args
    try:
        df = parse_argo_txt_file(input_path)
        if df is not None:
            output_path.parent.mkdir(parents=True, exist_ok=True)
            df.to_parquet(output_path, index=False)
    except Exception as e:
        print(f"\n[Error] 在工作进程中处理文件 {input_path.name} 时出错: {e}")


def process_argo_txt_to_yearly_parquet_dask(
    origin_dir: Path,
    temp_dir: Path,
    final_dir: Path,
    cleanup_temp_dir: bool = True
):
    """
    (主流程函数) 使用Dask统一并行处理Argo txt文件，高效处理大规模数据。

    功能:
        这是一个完整的、适用于超大数据集的ETL流程，完全由Dask驱动：
        1. (并行) [Dask] 将所有原始 .txt 文件转换为临时的 .parquet 文件。
        2. (并行) [Dask] 对所有临时文件进行并行地合并、清洗、排序，并输出为按年份分区的临时目录。
        3. (串行) 将分区目录合并为最终的单个年度文件。
        此版本利用Dask的HPC环境自适应能力，避免内存瓶颈，全程并行执行。
        最后，根据参数选择是否清理所有临时文件。

    参数:
        origin_dir (Path): 存放原始Argo .txt文件的目录。
        temp_dir (Path): 用于存放所有中间产物（初始Parquet、映射表、分区数据）的临时目录。
        final_dir (Path): 用于保存最终年份.parquet文件的目录。
        cleanup_temp_dir (bool, optional): 是否在任务结束后删除临时目录。默认为 True。
    """
    start_total_time = tm.time()
    
    # --- 准备工作：初始化Dask客户端并创建目录 ---
    client = Client()
    print(f"[*] Dask客户端已启动，仪表盘链接: {client.dashboard_link}")
    
    scheduler_info = client.scheduler_info()
    workers_count = len(scheduler_info['workers'])
    total_threads = sum(worker['nthreads'] for worker in scheduler_info['workers'].values())
    print(f"[*] Dask自动检测到 {workers_count} 个 worker，共 {total_threads} 个工作核心。")
    
    print("[*] 准备工作环境...")
    if temp_dir.exists():
        print(f"  - 发现旧的临时目录，正在清理: {temp_dir}")
        shutil.rmtree(temp_dir)
    temp_dir.mkdir(parents=True, exist_ok=True)
    final_dir.mkdir(parents=True, exist_ok=True)
    
    # --- 阶段一: 使用Dask并行解析 -> 中间文件 ---
    print("\n--- 阶段一: 开始并行解析原始 .txt 文件 ---")
    txt_files = list(origin_dir.glob('*.txt'))
    if not txt_files:
        print("[Error] 输入目录中未找到任何 .txt 文件。程序终止。")
        client.close()
        return
        
    tasks = [(p, temp_dir / f"{p.stem}.parquet") for p in txt_files]
    print(f"[*] 找到 {len(txt_files)} 个文件，将任务分发给Dask处理...")
    
    start_map_time = tm.time()
    futures = client.map(worker_process_file, tasks)
    client.gather(futures)
    print(f"--- 阶段一完成 --- (耗时: {tm.time() - start_map_time:.2f} 秒)")
    
    # --- 阶段二: Dask并行ETL并输出分区文件 ---
    print("\n--- 阶段二: Dask开始并行处理中间文件 ---")
    start_reduce_time = tm.time()

    print(f"[*] Dask正在逻辑读取 {temp_dir} 中的Parquet文件...")
    ddf = dd.read_parquet(temp_dir / '*.parquet')
    record_count = len(ddf)
    print(f"[*] 数据逻辑加载完成，总计 {record_count:,} 条记录。")

    print("[*] 正在重命名列...")
    rename_map = {'WMO': 'Platform_number', 'Cycle': 'Profile_number', 'Lon': 'Longitude', 'Lat': 'Latitude'}
    ddf = ddf.rename(columns=rename_map)

    print("[*] 正在并行生成全局唯一的 Profile_number...")
    unique_profiles = ddf[['Platform_number', 'Profile_number']].drop_duplicates().compute()
    unique_profiles_sorted = unique_profiles.sort_values(by=['Platform_number', 'Profile_number'])
    unique_profiles_sorted['Global_Profile_Number'] = range(1, len(unique_profiles_sorted) + 1)
    mapping_file_path = temp_dir / 'profile_mapping.parquet'
    unique_profiles_sorted.to_parquet(mapping_file_path, index=False)
    mapping_ddf = dd.read_parquet(mapping_file_path)
    ddf = ddf.merge(mapping_ddf, on=['Platform_number', 'Profile_number'], how='left')
    ddf = ddf.drop(columns=['Profile_number']).rename(columns={'Global_Profile_Number': 'Profile_number'})

    print("[*] 正在优化并转换数据类型...")
    float64_cols = list(ddf.select_dtypes(include='float64').columns)
    if float64_cols:
        for col in float64_cols:
            ddf[col] = ddf[col].astype('float32')
    int_conversion_map = {'Year': 'Int16', 'Month': 'Int8', 'Day': 'Int8', 'Platform_number': 'Int64', 'Profile_number': 'Int64'}
    for col, dtype in int_conversion_map.items():
        if col in ddf.columns:
            ddf[col] = ddf[col].astype(dtype)
    
    print("[*] Dask开始执行所有计算并按年份分区写入临时目录...")
    ddf = ddf.dropna(subset=['Year'])
    partitioned_output_dir = temp_dir / 'dask_partitioned_output'
    ddf.to_parquet(partitioned_output_dir, write_index=False, partition_on=['Year'])
    
    end_reduce_time = tm.time()
    print(f"--- 阶段二完成 --- (耗时: {end_reduce_time - start_reduce_time:.2f} 秒)")
    
    # --- 阶段三: 将分区文件合并为最终产物 ---
    print("\n--- 阶段三: 开始合并分区文件并重排定序 ---")
    start_consolidation_time = tm.time()
    
    year_dirs = sorted([d for d in partitioned_output_dir.iterdir() if d.is_dir() and d.name.startswith('Year=')])
    for year_dir in year_dirs:
        year_str = year_dir.name.split('=')[1]
        year_int = int(year_str)
        final_single_file_path = final_dir / f"Argo{year_int}.parquet"
        print(f"  -> 合并年份 {year_int} -> {final_single_file_path}")
        
        year_df = pd.read_parquet(year_dir)
        year_df['Year'] = year_int
        year_df.sort_values(by=['Year', 'Month', 'Day', 'Platform_number'], inplace=True)
        
        start_cols = ['Year', 'Month', 'Day']
        end_cols = ['Profile_number', 'Platform_number']
        middle_cols = [col for col in year_df.columns if col not in start_cols + end_cols]
        final_column_order = start_cols + middle_cols + end_cols
        year_df = year_df[final_column_order]
        
        year_df.to_parquet(final_single_file_path, index=False)
        
    end_consolidation_time = tm.time()
    print(f"--- 阶段三完成 --- (耗时: {end_consolidation_time - start_consolidation_time:.2f} 秒)")

    # --- 清理工作 ---
    if cleanup_temp_dir:
        print("\n--- 清理阶段: 删除临时文件 ---")
        try:
            shutil.rmtree(temp_dir)
            print(f"[*] 临时目录 {temp_dir} 已成功删除。")
        except Exception as e:
            print(f"[Error] 删除临时目录失败: {e}")
    else:
        print("\n--- 清理阶段: 已跳过 ---")
        print(f"[*] 临时文件已保留在: {temp_dir}")

    client.close()
    end_total_time = tm.time()
    print("\n==================================================")
    print(f"[Success] 所有任务完成！总耗时: {(end_total_time - start_total_time)/60:.2f} 分钟。")
    print("==================================================")

def load_argo_data(year: int, data_dir: str = './Argo_data', 
    variable_selection: dict | None = None) -> pd.DataFrame:
    """
    加载指定年份的 Argo Parquet 数据文件，并进行列名规范化和变量选择。

    功能:
        1. 自动将旧版列名 (如 'Depth_m') 转换为新版标准名 ('Depth')。
        2. 通过 `variable_selection` 参数，允许用户灵活选择最终输出的标准变量
           (Temperature, DO, Salinity) 分别来源于文件中的哪一列数据。

    参数:
        year (int): 
            需要加载的数据年份, 例如 2014。
        data_dir (str, optional): 
            存放 Argo Parquet 文件的目录。默认为 './Argo_data'。
        variable_selection (dict | None, optional): 
            一个字典，用于覆盖默认的变量来源。未指定的键将使用默认值。
            例如: {'Salinity': 'PSAL_WOA'} 只会更改盐度来源。
            默认来源为: {'Temperature': 'Temp_Adjusted', 'DO': 'DOXY_Adjusted', 'Salinity': 'PSAL_Adjusted'}

    返回:
        pd.DataFrame: 一个包含处理后 Argo 数据的 pandas DataFrame，其列名和数据源
                      均已根据参数进行了标准化。
    """
    # --- 1. 定义并合并变量选择 ---
    # 定义默认选择
    default_selection = {
        'Temperature': 'Temp_Adjusted',
        'Temperature_Flag': 'Temp_Adjusted_Flag',
        'DO': 'DOXY_Adjusted',
        'DO_Flag': 'DOXY_Adjusted_Flag',
        'Salinity': 'PSAL_Adjusted',
        'Salinity_Flag': 'PSAL_Adjusted_Flag'
    }
    # 如果用户提供了自定义选择，则用它来更新（覆盖）默认值
    if variable_selection:
        default_selection.update(variable_selection)
    
    # 最终生效的选择
    final_selection = default_selection
    
    # --- 2. 构建路径并加载文件 ---
    file_path = Path(data_dir) / f'Argo{year}.parquet'
    print(f"Attempting to load Argo data from: {file_path}")
    if not file_path.exists():
        raise FileNotFoundError(f"Error: The file '{file_path}' was not found.")
    try:
        argo_df = pd.read_parquet(file_path)
    except Exception as e:
        print(f"An error occurred while reading the Parquet file: {e}")
        raise

    # --- 3. 列名规范化：将所有已知的旧列名统一为新版中对应的名称 ---
    normalization_map = {
        'Depth_m': 'Depth',
        'Temperature_degC': 'Temp_Adjusted',
        'DO_mol_kg': 'DOXY_Adjusted',
        'Salinity_psu': 'PSAL_Adjusted'
    }
    argo_df.rename(columns=normalization_map, inplace=True)

    # --- 4. 变量选择与最终DataFrame构建 ---
    base_columns = [
        'Year', 'Month', 'Day', 'Longitude', 'Latitude', 'Depth',
        'Profile_number', 'Platform_number'
    ]
    existing_base_columns = [col for col in base_columns if col in argo_df.columns]
    final_df = argo_df[existing_base_columns].copy()

    for standard_name, source_col in final_selection.items():
        if source_col in argo_df.columns:
            final_df[standard_name] = argo_df[source_col]
        else:
            print(f"Warning: Column '{source_col}' not found in {file_path}. Creating empty column '{standard_name}'.")
            final_df[standard_name] = pd.NA

    print("Argo data loaded and processed successfully.")
    return final_df

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

def is_point_in_contour(row: pd.Series) -> bool:
    """
    检查一个数据行中的点坐标是否在其涡旋轮廓坐标内部。

    功能:
        这是一个辅助函数，通常与DataFrame的.apply()方法配合使用。
        它从一行数据中提取'Longitude', 'Latitude', 'contour_lon', 'contour_lat'，
        并判断点是否在多边形内部。

    参数:
        row (pd.Series): 
            一个Pandas DataFrame的数据行，必须包含经纬度和轮廓坐标。

    返回:
        bool: 如果点在轮廓内部，则返回True，否则返回False。
    """
    try:
        # 从行数据中提取轮廓坐标和点坐标
        contour_coords = zip(row['contour_lon'], row['contour_lat'])
        point_coord = Point(row['Longitude'], row['Latitude'])
        
        # 创建多边形对象
        polygon = Polygon(contour_coords)
        
        # 执行包含判断
        return polygon.contains(point_coord)
    except Exception:
        # 如果轮廓数据格式错误或不完整，安全地返回False
        return False

def filtered_float_data(
    DS: list, 
    no: int,
    argo_data_dir: str | Path = './Argo_data',
    circle_enlargement_factor: float = 1.2
) -> pd.DataFrame:
    """
    根据涡旋轨迹，动态加载并筛选匹配的Argo浮标剖面数据。

    功能:
        1. 分析涡旋轨迹覆盖的年份，并自动加载对应年份的Argo数据。
        2. 使用向量化的方式高效匹配在同一天出现在涡旋区域内的Argo浮标。
        3. 筛选标准为：浮标位置处于涡旋的有效轮廓内，或处于扩大一定倍数后的有效半径内。

    参数:
        DS (list): 
            包含所有涡旋轨迹信息的数据集。
        no (int): 
            需要筛选的涡旋的唯一编号。
        argo_data_dir (str | Path, optional): 
            存放Argo Parquet文件的目录。默认为 './Argo_data'。
        circle_enlargement_factor (float, optional): 
            筛选过程中涡旋半径的放大倍数。默认为 1.2。

    返回:
        pd.DataFrame: 一个包含所有匹配的Argo剖面完整数据的DataFrame（所有深度层级）。
                      如果无匹配数据，则返回一个空的DataFrame。
    """
    # --- 1. 准备涡旋数据 ---
    # print(f"[*] Preparing track data for eddy ID {no}...")
    wanted_track = find_track(DS, no)
    if not wanted_track:
        print(f"  - Track for eddy {no} not found, returning empty result.")
        return pd.DataFrame()

    # 将涡旋轨迹列表转换为DataFrame，方便后续处理
    track_df = pd.DataFrame(
        wanted_track,
        columns=['index_org', 'time', 'center_lon', 'center_lat', 'max_lon', 'max_lat', 
                 'contour_lon', 'contour_lat', 'radius', 'speed_contour_lon', 'speed_contour_lat']
    )
    # 使用convert_date函数转换日期
    track_df['date'] = convert_date(track_df['time'])

    # --- 2. 动态加载所需的Argo数据 ---
    min_year, max_year = track_df['date'].min().year, track_df['date'].max().year
    print(f"[*] Eddy track spans years: {min_year} - {max_year}. Loading corresponding Argo data...")
    
    all_argo_data = []
    for year in range(min_year, max_year + 1):
        try:
            yearly_argo_data = load_argo_data(year, data_dir=argo_data_dir)
            if not yearly_argo_data.empty:
                all_argo_data.append(yearly_argo_data)
        except FileNotFoundError:
            print(f"  - Warning: Argo data file for year {year} not found, skipping.")

    if not all_argo_data:
        print(f"  - No Argo data was loaded, returning empty result.")
        return pd.DataFrame()
        
    argo_data = pd.concat(all_argo_data, ignore_index=True)
    # print(f"[*] All relevant yearly Argo data loaded, total records: {len(argo_data):,}.")

    # --- 3. 准备Argo数据用于合并 ---
    # 预处理：删除没有有效坐标的记录，并创建日期列
    argo_data.dropna(subset=['Longitude', 'Latitude'], inplace=True)
    argo_data['date'] = pd.to_datetime(argo_data[['Year', 'Month', 'Day']])
    
    # 每个剖面只需要一个位置点来做匹配，因此我们按Profile_number去重，只保留第一个出现的点
    argo_positions = argo_data.drop_duplicates(subset='Profile_number', keep='first').copy()
    
    # --- 4. 核心步骤：将Argo与涡旋数据按日期合并 ---
    # pd.merge会高效地找出在同一天既有涡旋轨迹、又有Argo浮标的记录
    merged_df = pd.merge(argo_positions, track_df, on='date')
    if merged_df.empty:
        print(f"  - No Argo floats found on the same day as the eddy track.")
        return pd.DataFrame()
    # print(f"[*] Found {len(merged_df)} potential 'float-eddy' pairs based on date matching.")

    # --- 5. 向量化地理筛选 ---
    # print(f"[*] Performing geographic location matching...")
    # 5.1 检查是否在多边形内部
    inside_poly_mask = merged_df.apply(is_point_in_contour, axis=1)

    # 5.2 检查是否在扩大后的圆内（完全向量化）
    argo_coords = merged_df[['Longitude', 'Latitude']].values
    eddy_centers = merged_df[['center_lon', 'center_lat']].values
    distances = np.linalg.norm(argo_coords - eddy_centers, axis=1) # 批量计算所有点对的距离
    radii_deg = (merged_df['radius'].values / 111320) * circle_enlargement_factor
    inside_circle_mask = distances <= radii_deg
    
    # 5.3 合并两种筛选条件
    final_mask = inside_poly_mask | inside_circle_mask
    
    # 筛选出真正匹配的记录
    matched_profiles = merged_df[final_mask]
    
    if matched_profiles.empty:
        print(f"  - After geographic filtering, no matching Argo floats were found.")
        return pd.DataFrame()

    # --- 6. 根据筛选结果获取完整的剖面数据 ---
    # 从原始的、包含所有深度数据的argo_data中，筛选出所有匹配上的Profile_number
    matching_profile_numbers = matched_profiles['Profile_number'].unique()
    final_argo_data = argo_data[argo_data['Profile_number'].isin(matching_profile_numbers)].copy()
    
    # 在返回结果前，删除临时的'date'辅助列
    if 'date' in final_argo_data.columns:
        final_argo_data = final_argo_data.drop(columns=['date'])
    
    print(f"[*] Filtering complete! Found {len(matching_profile_numbers)} matching Argo profiles, with a total of {len(final_argo_data):,} data records.")
    
    return final_argo_data

def plot_track(
    DS: list, 
    no: int,
    save_fig: bool = False,
    show_fig: bool = True,
    plot_radius: bool = False,
    connection_threshold_days: int = 5,
    do_threshold: float = 50.0,
    salinity_threshold: float = 0.0,
    temperature_threshold: float = 0.0,
    depth_interval: float = 100.0,
    depth_merge_tolerance: float = 10.0,
    duplicate_depth_strategy: str = 'best_qc',
    anomaly_min_depth: float | None = 300.0,
    plot_unrelated_argo: bool = True,
    fix_delta_do_colorbar: bool = True,
    delta_do_cbar_min: float = 50.0,
    delta_do_cbar_max: float = 100.0,
    delta_do_cbar_ticks: list | None = None,
    min_anomaly_count: int = 0
):
    """
    绘制指定编号涡旋的详细轨迹，并智能高亮显示与 Argo 剖面的 ΔDO 异常交互情况。

    功能:
        1. 自动加载并筛选与指定涡旋匹配的Argo数据。
        2. 绘制涡旋的完整轨迹（虚线）。
        3. 智能高亮Argo浮标存在的时期：当浮标存在日的间隔小于阈值时，
           会将这段完整的涡旋轨迹绘制为连续实线（孤立的单日则标记为点）。
        4. 使用 calculate_delta_do 识别 ΔDO 异常（取每个剖面最大 ΔDO 一条，支持按最小深度过滤），并按 ΔDO 着色。
        5. 可选地绘制涡旋在交互日的有效半径和轮廓。

    参数:
        DS (list): 
            包含所有涡旋轨迹信息的数据集。
        no (int): 
            需要绘制的涡旋的唯一编号。
        save_fig (bool, optional): 
            是否将生成的图像保存到文件。默认为 False。
        show_fig (bool, optional): 
            是否在交互式环境中显示生成的图像。默认为 True。
        plot_radius (bool, optional): 
            是否以圆的形式绘制涡旋在交互日的有效半径。默认为 False。
        connection_threshold_days (int, optional):
            连接 Argo 交互点的最大天数阈值。默认为 5 天。
        do_threshold / salinity_threshold / temperature_threshold / depth_interval / depth_merge_tolerance / duplicate_depth_strategy: 
            传递给 calculate_delta_do 的参数。
        anomaly_min_depth (float | None): 
            ΔDO 异常最小深度限制；None 表示不限制。
        plot_unrelated_argo (bool): 
            是否绘制被 ΔDO 筛选掉的所有匹配 Argo 剖面基准位置（空心灰圈）。
        fix_delta_do_colorbar (bool): 
            是否固定 ΔDO 色标范围。
        delta_do_cbar_min / delta_do_cbar_max / delta_do_cbar_ticks: 
            色标范围与刻度设置。
        min_anomaly_count (int): 
            若 >0：要求 ΔDO 异常数量 ≥ 该值才绘图；=0 表示不做数量阈值过滤（默认 0）。
    """
    # --- 1. 准备涡旋和Argo数据 ---
    print(f"[*] Preparing data for eddy ID {no}...")
    
    # 获取涡旋轨迹并转换为DataFrame
    wanted_track = find_track(DS, no)
    if not wanted_track:
        print(f"  - Error: Track for eddy {no} not found.")
        return
        
    track_df = pd.DataFrame(
        wanted_track,
        columns=['index_org', 'time', 'center_lon', 'center_lat', 'max_lon', 'max_lat', 
                 'contour_lon', 'contour_lat', 'radius', 'speed_contour_lon', 'speed_contour_lat']
    )
    track_df['date'] = convert_date(track_df['time'])
    num = track_df['index_org'].iloc[0]

    # 获取数据集名称
    ds_name = "UNKNOWN"
    for name, var in inspect.currentframe().f_back.f_locals.items():
        if var is DS:
            ds_name = name.upper()
            break

    # 调用筛选函数，获取所有匹配的 Argo 数据（包含所有深度）
    argo_data_filtered = filtered_float_data(DS, no)

    # 预筛选：若匹配到的剖面数不足 min_anomaly_count，直接跳过
    if argo_data_filtered.empty:
        print(f"  - Skip plotting {ds_name}{no}: no matched Argo profiles.")
        return
    matched_profile_count = argo_data_filtered['Profile_number'].nunique()
    if min_anomaly_count > 0 and matched_profile_count < min_anomaly_count:
        print(f"\033[31m  - Skip plotting {ds_name}{no}: matched profiles ({matched_profile_count}) < min_anomaly_count ({min_anomaly_count}).\033[0m")
        return

    # 使用 ΔDO 异常检测
    anomalies = pd.DataFrame()
    base_argo_positions = pd.DataFrame()
    if not argo_data_filtered.empty:
        # 基准剖面位置（去除深度重复）
        base_argo_positions = (
            argo_data_filtered.sort_values(['Profile_number','Depth'])
            .groupby('Profile_number', as_index=False)
            .first()[['Profile_number','Longitude','Latitude','Year','Month','Day']]
        )
        anomalies = calculate_delta_do(
            argo_data_filtered,
            depth_interval=depth_interval,
            do_threshold=do_threshold,
            salinity_threshold=salinity_threshold,
            temperature_threshold=temperature_threshold,
            depth_merge_tolerance=depth_merge_tolerance,
            duplicate_depth_strategy=duplicate_depth_strategy,
            remove_outliers=True,
            verbose=False
        )
        if not anomalies.empty and anomaly_min_depth is not None:
            anomalies = anomalies[anomalies['depth'] >= anomaly_min_depth]
        if not anomalies.empty:
            # 每个剖面取最大 ΔDO
            anomalies = (
                anomalies.sort_values('delta_do', ascending=False)
                .drop_duplicates(subset='Profile_number', keep='first')
            )

    # 若异常数量不足阈值，直接跳过
    anomaly_count = 0 if anomalies.empty else len(anomalies)
    if min_anomaly_count > 0 and anomaly_count < min_anomaly_count:
        print(f"\033[31m  - Skip plotting {ds_name}{no}: anomalies ({anomaly_count}) < min_anomaly_count ({min_anomaly_count}).\033[0m")
        return

    # --- 2. 准备绘图 ---
    prop_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    eddy_color = prop_colors[1] if 'AC' in ds_name else prop_colors[0]
    
    # 加载地理底图
    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))

    # 初始化画布与坐标轴
    fig, ax = plt.subplots(figsize=(30, 20))
    ax.set_title(f'Track Analysis for Eddy {ds_name}{num}', fontsize=20)
    ax.set_xlabel('Longitude', fontsize=20); ax.set_ylabel('Latitude', fontsize=20)
    world.plot(color='lightgrey', edgecolor='white', ax=ax)

    # --- 3. 绘制涡旋轨迹 (背景和高亮) ---
    # 绘制完整的背景轨迹 (虚线)
    ax.plot(track_df['center_lon'], track_df['center_lat'], color=eddy_color, linestyle='--', alpha=0.5, label='Full Track Path')
    
    # 找出有Argo交互的日期
    interaction_track_df = pd.DataFrame(columns=track_df.columns)
    if not anomalies.empty:
        interaction_dates = pd.to_datetime(anomalies[['Year', 'Month', 'Day']]).unique()
        interaction_track_df = track_df[track_df['date'].isin(interaction_dates)].copy()
        
        # 识别并分别绘制不连续的实线段
        if not interaction_track_df.empty:
            date_diffs = interaction_track_df['date'].diff()
            break_points = date_diffs > pd.Timedelta(f'{connection_threshold_days} days')
            segment_ids = break_points.cumsum()
            
            labeled_interaction_path = False
            labeled_interaction_point = False
            # 按段落ID分组并分别绘图
            for _, segment_df in interaction_track_df.groupby(segment_ids):
                start_segment_date = segment_df['date'].min()
                end_segment_date = segment_df['date'].max()
                full_path_segment = track_df[
                    (track_df['date'] >= start_segment_date) & 
                    (track_df['date'] <= end_segment_date)
                ]

                # 如果段内只有一个交互点，画成散点
                if len(segment_df) == 1:
                    label = 'Interaction Point' if not labeled_interaction_point else None
                    ax.scatter(full_path_segment['center_lon'], full_path_segment['center_lat'], 
                               facecolors='none', edgecolors=eddy_color, 
                                linewidths=1.5, s=100, zorder=5, label=label)
                    labeled_interaction_point = True
                # 否则，将这段时间的完整轨迹画成实线
                else:
                    label = 'Interaction Path' if not labeled_interaction_path else None
                    ax.plot(full_path_segment['center_lon'], full_path_segment['center_lat'], 
                            color=eddy_color, linestyle='-', linewidth=2.5, label=label)
                    labeled_interaction_path = True
    
    # 标记轨迹的起点和终点
    ax.plot(track_df['center_lon'].iloc[0], track_df['center_lat'].iloc[0], marker='o', color=eddy_color, markersize=10, label='Start')
    ax.plot(track_df['center_lon'].iloc[-1], track_df['center_lat'].iloc[-1], marker='x', color=eddy_color, markersize=12, mew=2.5, label='End')

    # --- 4. 绘制交互日的轮廓、半径和日期 ---
    labeled_contour = False
    labeled_radius = False
    for _, eddy_day in interaction_track_df.iterrows():
        # 绘制轮廓线
        if not labeled_contour:
            ax.plot(eddy_day['contour_lon'], eddy_day['contour_lat'], color=eddy_color, linewidth=1, alpha=0.6, label='Effective Contour')
            labeled_contour = True
        else:
            ax.plot(eddy_day['contour_lon'], eddy_day['contour_lat'], color=eddy_color, linewidth=1, alpha=0.6)
        
        # 绘制有效半径
        if plot_radius:
            radius_color = 'r' if 'AC' in ds_name else 'purple'
            circle_label = 'Effective Radius' if not labeled_radius else None
            circle = plt.Circle((eddy_day['center_lon'], eddy_day['center_lat']), eddy_day['radius'] / 111320.0,
                                color=radius_color, fill=False, linestyle='--', alpha=0.4, linewidth=1.5, label=circle_label)
            ax.add_patch(circle)
            labeled_radius = True
            
        # 标记交互的起始和结束日期
        if eddy_day['date'] == interaction_track_df['date'].min() or eddy_day['date'] == interaction_track_df['date'].max():
             ax.text(eddy_day['center_lon'], eddy_day['center_lat'] + 0.1, eddy_day['date'].strftime('%Y-%m-%d'), 
                     fontsize=16, color='black', ha='center', zorder=11)

    # --- 5. 绘制 ΔDO 异常与基准剖面 ---
    if plot_unrelated_argo and not base_argo_positions.empty:
        ax.scatter(
            base_argo_positions['Longitude'], base_argo_positions['Latitude'],
            facecolors='none', edgecolors='gray', linewidths=0.8, s=50,
            label='All Matched Argo Profiles', zorder=5
        )

    if not anomalies.empty:
        scatter_kwargs = {}
        if fix_delta_do_colorbar:
            scatter_kwargs.update(dict(vmin=delta_do_cbar_min, vmax=delta_do_cbar_max))
        sc = ax.scatter(
            anomalies['Longitude'], anomalies['Latitude'],
            c=anomalies['delta_do'], cmap='Reds', s=90,
            edgecolors='black', linewidths=0.6,
            label=f'ΔDO ≥ {do_threshold} μmol kg⁻¹' + (f' @ depth ≥ {anomaly_min_depth} m' if anomaly_min_depth is not None else ''),
            zorder=10,
            **scatter_kwargs
        )
        cbar = plt.colorbar(sc, ax=ax, orientation='horizontal', fraction=0.046, pad=0.08)
        cbar.set_label('ΔDO / μmol·kg⁻¹', fontsize=18)
        cbar.ax.tick_params(labelsize=14)
        if fix_delta_do_colorbar:
            if delta_do_cbar_ticks is not None:
                cbar.set_ticks(delta_do_cbar_ticks)
            else:
                rng = delta_do_cbar_max - delta_do_cbar_min
                if rng > 30:
                    mid = (delta_do_cbar_min + delta_do_cbar_max) / 2
                    cbar.set_ticks([delta_do_cbar_min, mid, delta_do_cbar_max])
                else:
                    cbar.set_ticks([delta_do_cbar_min, delta_do_cbar_max])

    # --- 6. 最终化绘图设置 ---
    # 设定边界时排除META中错误的contour数据
    valid_contours_lon = [lon for day_lons in track_df['contour_lon'] for lon in day_lons if lon != 180.0]
    valid_contours_lat = [lat for day_lats in track_df['contour_lat'] for lat in day_lats if lat != 0.0]
    if valid_contours_lon and valid_contours_lat:
        ax.set_xlim(min(valid_contours_lon) - 0.5, max(valid_contours_lon) + 0.5)
        ax.set_ylim(min(valid_contours_lat) - 0.5, max(valid_contours_lat) + 0.5)
        
    ax.set_aspect('equal')
    ax.legend(fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=16)

    # --- 7. 输出控制 ---
    if save_fig:
        output_dir = Path("plot_track_analysis")
        output_dir.mkdir(exist_ok=True, parents=True)
        base_filename = f"Track_Analysis_{ds_name}{num}.png"
        save_path = output_dir / base_filename
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {save_path}")
    
    if show_fig:
        plt.show()
    
    plt.close(fig)
    
def convert_date(days_since_1950: pd.Series) -> pd.Series:
    """
    将以"自1950-01-01以来的天数"表示的数值转换为标准的datetime对象。

    功能:
        采用现代Pandas方法，使用pd.to_timedelta来处理日期运算，
        以避免版本更新带来的TypeError。

    参数:
        days_since_1950 (pd.Series): 包含天数数值的Pandas Series。

    返回:
        pd.Series: 转换后的datetime对象组成的Pandas Series。
    """
    # 定义基准日期，使用Pandas的Timestamp对象更佳
    t0 = pd.Timestamp('1950-01-01')
    
    # 核心步骤：将天数数值的Series转换为Timedelta Series
    # unit='D' 参数明确地告诉Pandas，这些数值的单位是“天”
    time_deltas = pd.to_timedelta(days_since_1950, unit='D')
    
    # 将基准日期与Timedelta Series相加，这是Pandas完全支持的操作
    return t0 + time_deltas

def plot_vertical(
    DS: list,
    no: int,
    show_fig: bool = False,
    save_fig: bool = False,
    color_mode: str = 'distance',
    variables: list = ['DO', 'Temp', 'Salinity'],
    show_colorbar: bool = False,
    remove_outliers: bool = True,
    aggregated: bool = False,
    argo_required: list | None = None,
    year_required: list | None = None,
    month_required: list | None = None,
    day_required: list | None = None
):
    '''
    根据涡旋轨迹与匹配到的 Argo 剖面，绘制变量-深度的垂直剖面。

    参数:
        DS (list): 涡旋轨迹数据集。
        no (int): 涡旋编号。
        show_fig (bool): 是否显示图片，默认 False。
        save_fig (bool): 是否保存图片，默认 False。
        color_mode (str): 颜色模式，'distance' 或 'time'，默认 'distance'。
        variables (list): 需要绘制的变量名称，默认 ['DO', 'Temp', 'Salinity']。
        show_colorbar (bool): 是否显示颜色条，默认 False。
        remove_outliers (bool): 是否执行 QC 过滤与规则法去极值，默认 True。
        aggregated (bool): 是否进行跨平台聚合绘制，默认 False。
        argo_required (list | None): 平台过滤；None 表示不过滤；传入平台编号列表时仅保留指定平台。
        year_required (list | None): 年份过滤；None 表示不过滤；传入年份列表时仅保留指定年份。
        month_required (list | None): 月份过滤；聚合模式 None 表示使用所有可用月份；逐平台模式 None 表示不过滤。
        day_required (list | None): 日期过滤（按日1-31）；None 表示不过滤；传入日期列表时仅保留指定日期。

    功能:
        - 为 variables 中的每个变量创建一个子图，按剖面绘制变量随深度变化的曲线。
        - 曲线颜色可根据与涡旋中心的相对距离（distance）或采样时间（time）变化。
        - 可选显示颜色条；支持图片保存与显示。
        - remove_outliers=True 时执行基础质量控制（QC 仅保留 {1,2,5,8}；DO<=1 置为 NaN）。
        - 可用 month_required 和 argo_required 对数据进行月份与平台筛选。

    模式差异:
        - aggregated=False（逐平台）：
            • 为每个浮标平台单独出图；每图包含 variables 中各变量的一个子图。
        - aggregated=True（聚合）：
            • 所有平台剖面在同一张图上聚合绘制（每个变量一个子图）。
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
        msg = "plot vertical profiles." if not aggregated else "to plot aggregated vertical profiles."
        print(f"No Argo data found for eddy {ds_names}{no} {msg}")
        return

    # 统一的月份和平台过滤（在两种模式下共享）
    # 先处理平台过滤（argo_required）
    if argo_required is not None:
        try:
            argo_required_set = set(int(x) for x in argo_required)
        except Exception:
            argo_required_set = set(argo_required)
        before_cnt = len(argo_data_filtered)
        argo_data_filtered = argo_data_filtered[argo_data_filtered['Platform_number'].isin(argo_required_set)]
        if argo_data_filtered.empty:
            print(f"No Argo data left after filtering by platforms {sorted(list(argo_required_set))}.")
            return

    # 处理年月日过滤（year/month/day）
    # 注意：aggregated=True 时，按原逻辑在未提供月份时默认使用所有可用月份；
    #       aggregated=False 时，仅在提供了 month_required 时进行过滤。
    #       year_required/day_required 在两种模式下仅当提供时进行过滤。
    
    # Year
    argo_data_filtered['Year'] = pd.to_numeric(argo_data_filtered['Year'], errors='coerce')
    argo_data_filtered.dropna(subset=['Year'], inplace=True)
    argo_data_filtered['Year'] = argo_data_filtered['Year'].astype(int)

    # Month
    argo_data_filtered['Month'] = pd.to_numeric(argo_data_filtered['Month'], errors='coerce')
    argo_data_filtered.dropna(subset=['Month'], inplace=True)
    argo_data_filtered['Month'] = argo_data_filtered['Month'].astype(int)

    # Day
    argo_data_filtered['Day'] = pd.to_numeric(argo_data_filtered['Day'], errors='coerce')
    argo_data_filtered.dropna(subset=['Day'], inplace=True)
    argo_data_filtered['Day'] = argo_data_filtered['Day'].astype(int)

    if aggregated:
        # 若未指定月份，采用所有可用月份（保持原 plot_vertical_monthly 行为）
        if not month_required:
            all_available_months = sorted(argo_data_filtered['Month'].unique().tolist())
            if not all_available_months:
                print(f"No valid months found in data for eddy {ds_names}{no}.")
                return
            print(f"No months specified. Defaulting to all available months: {all_available_months}")
            month_required = all_available_months
        # 根据年月日过滤
        if month_required:
            argo_data_filtered = argo_data_filtered[argo_data_filtered['Month'].isin(month_required)].copy()
        if year_required:
            argo_data_filtered = argo_data_filtered[argo_data_filtered['Year'].isin(year_required)].copy()
        if day_required:
            argo_data_filtered = argo_data_filtered[argo_data_filtered['Day'].isin(day_required)].copy()
        if argo_data_filtered.empty:
            print(f"No data found for eddy {ds_names}{no} after applying filters: year={year_required}, month={month_required}, day={day_required}.")
            return
    else:
        # 非聚合模式：仅当提供了对应过滤条件才过滤
        if year_required:
            argo_data_filtered = argo_data_filtered[argo_data_filtered['Year'].isin(year_required)].copy()
            if argo_data_filtered.empty:
                print(f"No data found for eddy {ds_names}{no} in years {year_required} (per-platform mode).")
                return
        if month_required:
            argo_data_filtered = argo_data_filtered[argo_data_filtered['Month'].isin(month_required)].copy()
            if argo_data_filtered.empty:
                print(f"No data found for eddy {ds_names}{no} in months {month_required} (per-platform mode).")
                return
        if day_required:
            argo_data_filtered = argo_data_filtered[argo_data_filtered['Day'].isin(day_required)].copy()
            if argo_data_filtered.empty:
                print(f"No data found for eddy {ds_names}{no} in days {day_required} (per-platform mode).")
                return

    # =========================
    # 分支一：非聚合（原 plot_vertical 行为，新增月份/平台筛选）
    # =========================
    if not aggregated:
        for platform_id_val, platform_data in argo_data_filtered.groupby("Platform_number"):
            profile_num_agg = platform_data['Profile_number'].agg(['min', 'max'])
            min_profile_num = profile_num_agg['min']
            max_profile_num = profile_num_agg['max']

            num_variables = len(variables)
            fig, axes = plt.subplots(1, num_variables, figsize=(10 * num_variables, 20), sharey=True)
            if num_variables == 1:
                axes = [axes]

            cmap = plt.cm.coolwarm
            profile_dates_for_title = []

            for i, var_name in enumerate(variables):
                ax = axes[i]
                original_variable_name = var_name
                db_variable_name = 'Temperature' if var_name == 'Temp' else var_name

                if db_variable_name not in platform_data.columns:
                    ax.text(0.5, 0.5, f"Variable '{db_variable_name}'\nnot found in data.",
                            ha='center', va='center', transform=ax.transAxes, fontsize=16)
                    ax.set_title(f"Variable: {original_variable_name}", fontsize=20)
                    continue

                for profile_num, rows in platform_data.groupby("Profile_number"):
                    if rows.empty:
                        continue

                    rows_to_plot = rows.dropna(subset=[db_variable_name, "Depth"])
                    if rows_to_plot.empty:
                        continue

                    if remove_outliers:
                        qc_column_name = f"{db_variable_name}_Flag"
                        if qc_column_name in rows_to_plot.columns:
                            good_qc_flags = ['1', '2', '5', '8', 1, 2, 5, 8]
                            bad_qc_mask = ~rows_to_plot[qc_column_name].isin(good_qc_flags)
                            rows_to_plot.loc[bad_qc_mask, db_variable_name] = np.nan
                        if db_variable_name == 'DO':
                            bad_value_mask = rows_to_plot[db_variable_name] <= 1.0
                            rows_to_plot.loc[bad_value_mask, db_variable_name] = np.nan

                    if rows_to_plot.empty:
                        continue

                    try:
                        current_profile_date = pd.Timestamp(year=int(rows.iloc[0]['Year']),
                                                            month=int(rows.iloc[0]['Month']),
                                                            day=int(rows.iloc[0]['Day']))
                        if i == 0:
                            profile_dates_for_title.append(current_profile_date)
                    except (ValueError, TypeError):
                        continue

                    color_value_normalized = 0.5

                    if color_mode == 'distance':
                        if 'Longitude' in rows.iloc[0] and 'Latitude' in rows.iloc[0] and wanted_track:
                            track_dates_converted = convert_date([t[1] for t in wanted_track])
                            idx_track_list = [j for j, td in enumerate(track_dates_converted) if hasattr(td, 'date') and td.date() == current_profile_date.date()]
                            if idx_track_list:
                                idx_track = idx_track_list[0]
                                center_lon, center_lat, radius = wanted_track[idx_track][2], wanted_track[idx_track][3], wanted_track[idx_track][8]
                                if radius > 1e-6:
                                    rel_x = (rows.iloc[0]['Longitude'] - center_lon) / (radius / 111320.0)
                                    rel_y = (rows.iloc[0]['Latitude'] - center_lat) / (radius / 111320.0)
                                    distance = np.sqrt(rel_x**2 + rel_y**2)
                                    color_value_normalized = 1.0 - np.clip(distance, 0.0, 1.0)
                    elif color_mode == 'time':
                        if max_profile_num > min_profile_num:
                            color_value_normalized = (profile_num - min_profile_num) / (max_profile_num - min_profile_num)
                        else:
                            color_value_normalized = 0.0

                    color = cmap(np.clip(color_value_normalized, 0.0, 1.0))
                    ax.plot(rows_to_plot[db_variable_name], rows_to_plot["Depth"], color=color, alpha=0.7)

                # 子图属性
                ax.set_ylim(-50, 2050)
                if db_variable_name == 'DO':
                    ax.set_xlim(10, 350)
                elif db_variable_name == 'Temperature':
                    ax.set_xlim(1, 32)
                elif db_variable_name == 'Salinity':
                    ax.set_xlim(32.5, 35.5)

                ax.set_xlabel(original_variable_name, fontsize=20)
                ax.tick_params(axis='x', labelsize=16)
                ax.grid(True)

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

            # 整体属性
            if not profile_dates_for_title:
                plt.close(fig)
                continue

            date_start_platform = min(profile_dates_for_title)
            date_end_platform = max(profile_dates_for_title)

            axes[0].set_ylabel("Depth/m", fontsize=20)
            axes[0].tick_params(axis='y', labelsize=16)
            axes[0].invert_yaxis()

            fig.suptitle(f"{ds_names}{no}, Platform: {int(platform_id_val)}, {date_start_platform.date()}~{date_end_platform.date()}", fontsize=24, y=0.95)
            plt.tight_layout(rect=[0, 0, 1, 0.93])

            if save_fig:
                output_dir = "plot_vertical_profiles"
                os.makedirs(output_dir, exist_ok=True)
                plt.savefig(os.path.join(output_dir, f"{ds_names}{no}_Platform_{int(platform_id_val)}.png"), dpi=300, bbox_inches='tight')

            if show_fig:
                plt.show()
            plt.close(fig)
        return

    # =========================
    # 分支二：聚合（等价于原 plot_vertical_monthly）
    # =========================
    # 收集各剖面的基本信息
    profiles_to_plot = []
    try:
        argo_data_filtered['date_ts'] = pd.to_datetime(argo_data_filtered[['Year', 'Month', 'Day']])
    except (ValueError, TypeError) as e:
        print(f"Could not create timestamps for all rows due to invalid date components: {e}.")
        argo_data_filtered = argo_data_filtered.dropna(subset=['Year', 'Month', 'Day'])
        argo_data_filtered['date_ts'] = pd.to_datetime(argo_data_filtered[['Year', 'Month', 'Day']])

    for _, profile_rows in argo_data_filtered.groupby(['Platform_number', 'Profile_number']):
        if not profile_rows.empty:
            profiles_to_plot.append({
                'rows': profile_rows,
                'date': profile_rows.iloc[0]['date_ts'],
                'lon': profile_rows.iloc[0]['Longitude'],
                'lat': profile_rows.iloc[0]['Latitude']
            })

    if not profiles_to_plot:
        print(f"No data found for eddy {ds_names}{no} in months {month_required}.")
        return

    all_dates = [p['date'] for p in profiles_to_plot]
    min_time_for_norm, max_time_for_norm = (min(all_dates), max(all_dates)) if all_dates else (None, None)

    num_variables = len(variables)
    fig, axes = plt.subplots(1, num_variables, figsize=(10 * num_variables, 20), sharey=True)
    if num_variables == 1:
        axes = [axes]

    cmap = plt.cm.coolwarm

    for i, var_name in enumerate(variables):
        ax = axes[i]
        original_variable_name = var_name
        db_variable_name = 'Temperature' if var_name == 'Temp' else var_name

        if db_variable_name not in argo_data_filtered.columns:
            ax.text(0.5, 0.5, f"Variable '{db_variable_name}'\nnot found in data.",
                    ha='center', va='center', transform=ax.transAxes, fontsize=16)
            ax.set_xlabel(original_variable_name, fontsize=20)
            ax.grid(True)
            continue

        for profile_info in profiles_to_plot:
            rows = profile_info['rows']
            rows_to_plot = rows.dropna(subset=[db_variable_name, "Depth"])
            if rows_to_plot.empty:
                continue

            if remove_outliers:
                qc_column_name = f"{db_variable_name}_Flag"
                if qc_column_name in rows_to_plot.columns:
                    good_qc_flags = ['1', '2', '5', '8', 1, 2, 5, 8]
                    bad_qc_mask = ~rows_to_plot[qc_column_name].isin(good_qc_flags)
                    rows_to_plot.loc[bad_qc_mask, db_variable_name] = np.nan
                if db_variable_name == 'DO':
                    bad_value_mask = rows_to_plot[db_variable_name] <= 1.0
                    rows_to_plot.loc[bad_value_mask, db_variable_name] = np.nan

            if rows_to_plot.empty:
                continue

            current_date = profile_info['date']
            color_value_normalized = 0.5

            if color_mode == 'distance':
                if wanted_track and 'lon' in profile_info and 'lat' in profile_info:
                    track_dates_converted = convert_date([t[1] for t in wanted_track])
                    idx_track_list = [j for j, td in enumerate(track_dates_converted) if hasattr(td, 'date') and td.date() == current_date.date()]
                    if idx_track_list:
                        idx_track = idx_track_list[0]
                        center_lon, center_lat, radius = wanted_track[idx_track][2], wanted_track[idx_track][3], wanted_track[idx_track][8]
                        if radius > 1e-6:
                            rel_x = (profile_info['lon'] - center_lon) / (radius / 111320.0)
                            rel_y = (profile_info['lat'] - center_lat) / (radius / 111320.0)
                            distance = np.sqrt(rel_x**2 + rel_y**2)
                            color_value_normalized = 1.0 - np.clip(distance, 0.0, 1.0)
            elif color_mode == 'time':
                if min_time_for_norm and max_time_for_norm and max_time_for_norm > min_time_for_norm:
                    total_delta = (max_time_for_norm - min_time_for_norm).total_seconds()
                    current_delta = (current_date - min_time_for_norm).total_seconds()
                    color_value_normalized = current_delta / total_delta if total_delta > 0 else 0.0
                else:
                    color_value_normalized = 0.0

            color = cmap(np.clip(color_value_normalized, 0.0, 1.0))
            ax.plot(rows_to_plot[db_variable_name], rows_to_plot["Depth"], color=color, alpha=0.7)

        # 子图属性
        ax.set_ylim(-50, 2050)
        if db_variable_name == 'DO':
            ax.set_xlim(10, 350)
        elif db_variable_name == 'Temperature':
            ax.set_xlim(1, 32)
        elif db_variable_name == 'Salinity':
            ax.set_xlim(32.5, 35.5)

        ax.set_xlabel(original_variable_name, fontsize=20)
        ax.tick_params(axis='x', labelsize=16)
        ax.grid(True)

        if show_colorbar:
            norm = Normalize(vmin=0, vmax=1)
            sm = ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = plt.colorbar(sm, ax=ax, orientation='vertical', fraction=0.046, pad=0.04)
            if color_mode == 'distance':
                cbar.set_ticks([0, 0.5, 1])
                cbar.set_ticklabels(['1.0', '0.5', '0.0'])
                cbar.set_label('Normalized Distance (0=center, 1=edge)', fontsize=14)
            elif color_mode == 'time' and min_time_for_norm and max_time_for_norm:
                cbar.set_label(f'Normalized Time\n({min_time_for_norm.strftime("%Y-%m-%d")} to {max_time_for_norm.strftime("%Y-%m-%d")})', fontsize=12)

    # 整体属性
    axes[0].set_ylabel("Depth/m", fontsize=20)
    axes[0].tick_params(axis='y', labelsize=16)
    axes[0].invert_yaxis()

    month_str = "All" if month_required and len(month_required) > 6 else (", ".join(map(str, month_required)) if month_required else "All")
    date_range_str = f"{min_time_for_norm.date()}~{max_time_for_norm.date()}" if all_dates else "No date range"
    fig.suptitle(f"{ds_names}{no}, Months: {month_str}, Data: {date_range_str}", fontsize=24, y=0.95)
    plt.tight_layout(rect=[0, 0, 1, 0.93])

    if save_fig:
        output_dir = "plot_vertical_monthly_aggregated"
        os.makedirs(output_dir, exist_ok=True)
        month_suffix = "all" if not month_required or (month_required and len(month_required) > 6) else "_".join(map(str, month_required))
        filename = f"{ds_names}{no}_months_{month_suffix}_aggregated.png"
        plt.savefig(os.path.join(output_dir, filename), dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {os.path.join(output_dir, filename)}")

    if show_fig:
        plt.show()
    plt.close(fig)
        
def plot_relative_position(
    DS: list,
    no: int,
    show_fig: bool = False,
    save_fig: bool = False,
    color_mode: str = 'distance',
    show_colorbar: bool = False,
    aggregated: bool = False,
    argo_required: list | None = None,
    year_required: list | None = None,
    month_required: list | None = None,
    day_required: list | None = None
):
    '''
    根据涡旋轨迹和浮标数据，绘制浮标在单位圆涡旋中的相对位置分布图。

    参数:
        DS (list): 涡旋轨迹数据集。
        no (int): 涡旋编号。
        show_fig (bool): 是否显示图片，默认 False。
        save_fig (bool): 是否保存图片，默认 False。
        color_mode (str): 'distance' 或 'time'，默认 'distance'。
        show_colorbar (bool): 是否显示颜色条，默认 False。
        aggregated (bool): 是否聚合所有平台于一图，默认 False。
        argo_required (list | None): 平台筛选；None 表示不过滤；传入平台编号列表时仅保留指定平台。
        year_required (list | None): 年份筛选；None 表示不过滤。
        month_required (list | None): 月份筛选；聚合模式 None 表示使用所有可用月份；逐平台模式 None 表示不过滤。
        day_required (list | None): 日期筛选（按日1-31）；None 表示不过滤。

        功能:
        - 计算剖面代表点在单位圆涡旋中的相对位置，并以散点标注（中心×，单位圆圈）。
        - 颜色模式: 'distance'（距离中心归一化，0=中心，1=边缘）或 'time'（按时间顺序归一化）。
        - 横纵轴刻度包含真实经纬度与相对坐标。可选显示颜色条；支持图片保存与显示。
        - 支持按照月份（month_required）与平台编号（argo_required）进行筛选。

    模式差异:
        - aggregated=False（逐平台）：
            • 对每个浮标平台分别出图，点内数字代表该平台内部的剖面时序，从1开始递增。
        - aggregated=True（聚合）：
            • 所有平台的代表点聚合到同一张图，点内数字代表相对于所选月份范围起始日的累积天数。例如，若数据从7月29日开始，则该日所有剖面的数字为29，7月30日为30，8月1日则为32。
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

    # 平台筛选（与 plot_vertical 一致）
    if argo_required is not None:
        try:
            argo_required_set = set(int(x) for x in argo_required)
        except Exception:
            argo_required_set = set(argo_required)
        argo_data_filtered = argo_data_filtered[argo_data_filtered['Platform_number'].isin(argo_required_set)]
        if argo_data_filtered.empty:
            print(f"No Argo data left after filtering by platforms {sorted(list(argo_required_set))}.")
            return

    # 年月日筛选（与 plot_vertical 一致）
    argo_data_filtered['Year'] = pd.to_numeric(argo_data_filtered['Year'], errors='coerce')
    argo_data_filtered.dropna(subset=['Year'], inplace=True)
    argo_data_filtered['Year'] = argo_data_filtered['Year'].astype(int)

    argo_data_filtered['Month'] = pd.to_numeric(argo_data_filtered['Month'], errors='coerce')
    argo_data_filtered.dropna(subset=['Month'], inplace=True)
    argo_data_filtered['Month'] = argo_data_filtered['Month'].astype(int)

    argo_data_filtered['Day'] = pd.to_numeric(argo_data_filtered['Day'], errors='coerce')
    argo_data_filtered.dropna(subset=['Day'], inplace=True)
    argo_data_filtered['Day'] = argo_data_filtered['Day'].astype(int)

    if aggregated:
        if not month_required:
            all_months = sorted(argo_data_filtered['Month'].unique().tolist())
            if not all_months:
                print(f"No valid months found in data for eddy {ds_names}{no}.")
                return
            print(f"No months specified for relative position plot. Defaulting to all available months: {all_months}")
            month_required = all_months
        # 按提供的年/月/日过滤（月份在聚合模式下为必有列表）
        if month_required:
            argo_data_filtered = argo_data_filtered[argo_data_filtered['Month'].isin(month_required)].copy()
        if year_required:
            argo_data_filtered = argo_data_filtered[argo_data_filtered['Year'].isin(year_required)].copy()
        if day_required:
            argo_data_filtered = argo_data_filtered[argo_data_filtered['Day'].isin(day_required)].copy()
        if argo_data_filtered.empty:
            print(f"No data found for eddy {ds_names}{no} after applying filters: year={year_required}, month={month_required}, day={day_required}.")
            return
    else:
        if year_required:
            argo_data_filtered = argo_data_filtered[argo_data_filtered['Year'].isin(year_required)].copy()
            if argo_data_filtered.empty:
                print(f"No data found for eddy {ds_names}{no} in years {year_required} (per-platform mode).")
                return
        if month_required:
            argo_data_filtered = argo_data_filtered[argo_data_filtered['Month'].isin(month_required)].copy()
            if argo_data_filtered.empty:
                print(f"No data found for eddy {ds_names}{no} in months {month_required} (per-platform mode).")
                return
        if day_required:
            argo_data_filtered = argo_data_filtered[argo_data_filtered['Day'].isin(day_required)].copy()
            if argo_data_filtered.empty:
                print(f"No data found for eddy {ds_names}{no} in days {day_required} (per-platform mode).")
                return

    # =========================
    # 分支一：逐平台模式（aggregated=False）
    # =========================
    if not aggregated:
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

            # 设置坐标轴刻度 (与聚合版一致)
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
        return

    # =========================
    # 分支二：聚合模式（aggregated=True）
    # =========================
    # 1) 收集聚合绘图所需点
    points_to_process = []  # {date, data_row}
    monthly_filtered_data = argo_data_filtered.copy()
    if month_required:
        monthly_filtered_data = monthly_filtered_data[monthly_filtered_data['Month'].isin(month_required)].copy()

    if not monthly_filtered_data.empty:
        # 获取每个剖面的第一行数据作为代表点
        profile_first_rows = monthly_filtered_data.groupby(["Platform_number", "Profile_number"]).first().reset_index()
        for _, p_row in profile_first_rows.iterrows():
            try:
                current_date_profile = pd.Timestamp(year=int(p_row['Year']), month=int(p_row['Month']), day=int(p_row['Day']))
            except (ValueError, TypeError):
                continue
            if not month_required or current_date_profile.month in month_required:
                points_to_process.append({'date': current_date_profile, 'data_row': p_row})

    if not points_to_process:
        print(f"No data found for eddy {ds_names}{no} in months {month_required}.")
        return

    # 2) 计算日期标签的参考起点（按所选月份的最小月的一号）
    min_plot_date_overall = min(p['date'] for p in points_to_process)
    ref_month = min(month_required) if month_required else min_plot_date_overall.month
    reference_start_date_for_labels = pd.Timestamp(year=min_plot_date_overall.year, month=ref_month, day=1)

    points_to_plot = []
    all_track_info_for_overall_mean = []
    all_profile_dates_for_title = []
    all_profile_timestamps_for_time_mode = []

    track_dates_converted = convert_date([t[1] for t in wanted_track]) if wanted_track else []

    for point_info in points_to_process:
        current_date = point_info['date']
        p_row = point_info['data_row']

        day_label = (current_date - reference_start_date_for_labels).days + 1

        center_lon, center_lat, radius = None, None, None
        if wanted_track:
            matches = [i for i, td in enumerate(track_dates_converted) if hasattr(td, 'date') and td.date() == current_date.date()]
            if matches:
                idx_track = matches[0]
                center_lon = wanted_track[idx_track][2]
                center_lat = wanted_track[idx_track][3]
                radius = wanted_track[idx_track][8]

        if center_lon is not None and radius is not None and radius > 1e-6:
            if 'Longitude' not in p_row or 'Latitude' not in p_row:
                continue

            rel_x = (p_row['Longitude'] - center_lon) / (radius / 111320.0)
            rel_y = (p_row['Latitude'] - center_lat) / (radius / 111320.0)

            points_to_plot.append({'rel_x': rel_x, 'rel_y': rel_y, 'date': current_date, 'day_label': day_label})
            all_track_info_for_overall_mean.append([center_lon, center_lat, radius])
            all_profile_dates_for_title.append(current_date)
            if color_mode == 'time':
                all_profile_timestamps_for_time_mode.append(current_date)

    if not points_to_plot:
        print(f"No valid points with track data found for eddy {ds_names}{no} in months {month_required}.")
        return

    # 3) 确定时间归一化范围（仅用于 color_mode='time'）
    min_time_for_norm, max_time_for_norm = (None, None)
    if color_mode == 'time' and all_profile_timestamps_for_time_mode:
        min_time_for_norm = min(all_profile_timestamps_for_time_mode)
        max_time_for_norm = max(all_profile_timestamps_for_time_mode)
        if min_time_for_norm == max_time_for_norm and len(all_profile_timestamps_for_time_mode) > 1:
            max_time_for_norm = min_time_for_norm + pd.Timedelta(days=1)

    date_start_overall = min(all_profile_dates_for_title) if all_profile_dates_for_title else None
    date_end_overall = max(all_profile_dates_for_title) if all_profile_dates_for_title else None

    # 4) 绘图
    fig, ax = plt.subplots(figsize=(30, 20))
    cmap = plt.cm.coolwarm

    for point in points_to_plot:
        rel_x, rel_y = point['rel_x'], point['rel_y']
        current_date, day_label = point['date'], point['day_label']

        color_value_normalized = 0.5
        if color_mode == 'distance':
            distance_from_center = np.sqrt(rel_x**2 + rel_y**2)
            color_value_normalized = 1.0 - np.clip(distance_from_center, 0.0, 1.0)
        elif color_mode == 'time':
            if min_time_for_norm and max_time_for_norm and min_time_for_norm < max_time_for_norm:
                total_delta = (max_time_for_norm - min_time_for_norm).total_seconds()
                current_delta = (current_date - min_time_for_norm).total_seconds()
                color_value_normalized = (current_delta / total_delta) if total_delta > 0 else 0.0
            elif min_time_for_norm and max_time_for_norm and min_time_for_norm == max_time_for_norm:
                color_value_normalized = 0.0

        color = cmap(np.clip(color_value_normalized, 0.0, 1.0))
        ax.scatter(rel_x, rel_y, color=color, s=300, zorder=5)
        ax.text(rel_x, rel_y, str(day_label), weight='bold', fontsize=9, color='black', ha='center', va='center', zorder=6)

    # 中心与单位圆
    ax.plot(0, 0, marker='x', color='black', markersize=16, markeredgewidth=3, label='Eddy Center (Relative)', zorder=3)
    circle = plt.Circle((0, 0), 1, color='gray', fill=False, linestyle='--', linewidth=2, label='Unit Eddy Boundary', zorder=2)
    ax.add_patch(circle)
    ax.set_aspect('equal')

    # 标题与坐标轴
    month_str = "All" if (month_required and len(month_required) > 6) else ", ".join(map(str, month_required)) if month_required else "All"
    title_str = f"{ds_names}{no}, Months: {month_str}, Relative Positions"
    if date_start_overall and date_end_overall:
        title_str += f"\nData: {date_start_overall.date()}~{date_end_overall.date()}, Total Points: {len(points_to_plot)}"
    ax.set_title(title_str, fontsize=20)
    ax.set_xlabel('Relative X (Eddy Radii)', fontsize=20)
    ax.set_ylabel('Relative Y (Eddy Radii)', fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=16)

    # 地理刻度（均值）
    if all_track_info_for_overall_mean:
        mean_center_lon = np.mean([info[0] for info in all_track_info_for_overall_mean])
        mean_center_lat = np.mean([info[1] for info in all_track_info_for_overall_mean])
        mean_radius = np.mean([info[2] for info in all_track_info_for_overall_mean])
        if not np.isnan(mean_center_lon) and not np.isnan(mean_center_lat) and not np.isnan(mean_radius) and mean_radius > 1e-6:
            mean_degrees = mean_radius / 111320.0
            tick_locs = [-1, -0.5, 0, 0.5, 1]
            x_tick_labels = [f"{(mean_center_lon + t * mean_degrees):.2f}°\n({t})" for t in tick_locs]
            y_tick_labels = [f"{(mean_center_lat + t * mean_degrees):.2f}°\n({t})" for t in tick_locs]
            ax.set_xticks(tick_locs)
            ax.set_xticklabels(x_tick_labels)
            ax.set_yticks(tick_locs)
            ax.set_yticklabels(y_tick_labels)

    ax.set_xlim([-1.25, 1.25])
    ax.set_ylim([-1.25, 1.25])

    # 颜色条
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
            if min_time_for_norm and max_time_for_norm:
                cbar.set_label(f'Normalized Time ({min_time_for_norm.strftime("%Y-%m-%d")} to {max_time_for_norm.strftime("%Y-%m-%d")})', fontsize=14)
            else:
                cbar.set_label('Normalized Time', fontsize=14)

    # 保存/显示
    if save_fig:
        output_dir = "plot_relative_position_monthly_aggregated"
        os.makedirs(output_dir, exist_ok=True)
        month_suffix = "all" if not month_required or (month_required and len(month_required) > 6) else "_".join(map(str, month_required))
        base_filename = f"{ds_names}{no}_RP_months_{month_suffix}_aggregated.png"
        plt.savefig(os.path.join(output_dir, base_filename), dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {os.path.join(output_dir, base_filename)}")

    if show_fig:
        plt.show()
    plt.close(fig)

def plot_vertical_monthly(
    DS: list,
    no: int,
    month_required: list = None,
    show_fig: bool = False,
    save_fig: bool = False,
    color_mode: str = 'distance',
    variables: list = ['DO', 'Temp', 'Salinity'],
    show_colorbar: bool = False,
    remove_outliers: bool = True,
    argo_required: list | None = None
):
    '''
    兼容封装：保持旧接口，内部调用统一的 plot_vertical。
    '''
    return plot_vertical(
        DS=DS,
        no=no,
        show_fig=show_fig,
        save_fig=save_fig,
        color_mode=color_mode,
        variables=variables,
        show_colorbar=show_colorbar,
        remove_outliers=remove_outliers,
        aggregated=True,
        month_required=month_required,
        argo_required=argo_required
    )
    
def plot_relative_position_monthly(
    DS: list,
    no: int,
    show_fig: bool = False,
    save_fig: bool = False,
    color_mode: str = 'distance',
    show_colorbar: bool = False,
    month_required: list | None = None,
    argo_required: list | None = None
):
    '''
    兼容封装：保持旧接口，内部调用统一的 plot_relative_position（aggregated=True）。
    '''
    return plot_relative_position(
        DS=DS,
        no=no,
        show_fig=show_fig,
        save_fig=save_fig,
        color_mode=color_mode,
        show_colorbar=show_colorbar,
        aggregated=True,
        month_required=month_required,
        argo_required=argo_required,
    )

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

class LineDrawer:
    """
    一个用于在matplotlib图像上通过点击两点来交互式绘制直线的辅助类。
    新增了撤回功能：按 'z' 键可以撤回上一个点击的点或已绘制的直线。
    """
    def __init__(self, ax):
        self.ax = ax
        self.points = []  # 用于存储用户点击的坐标
        self.lines = []   # 用于存储已绘制的直线对象
        self.markers = [] # 用于存储用户点击时产生的临时视觉标记
        # 将键盘事件与 onkey 方法连接
        self.cid_key = self.ax.figure.canvas.mpl_connect('key_press_event', self.onkey)
        
        print("\n--- 交互模式已激活 (稳定版) ---")
        print("请在图像上点击两个点来定义一条新的直线。")
        print("提示：按 'z' 键可以撤回上一个点或刚绘制的直线。")

    def onclick(self, event):
        """鼠标点击事件的处理函数。"""
        # 确保点击事件发生在指定的坐标轴内
        if event.inaxes != self.ax: return
        
        # 记录点击的坐标点
        self.points.append((event.xdata, event.ydata))
        # 在点击位置添加一个临时的"+"标记，给予用户视觉反馈
        marker = self.ax.plot(event.xdata, event.ydata, 'm+', markersize=12, markeredgewidth=2)
        self.markers.extend(marker)
        # 触发一次完整的重绘，确保标记点被永久画上
        self.ax.figure.canvas.draw()

        # 当记录的点达到两个时，开始绘制直线
        if len(self.points) == 2:
            self.draw_line()

    def onkey(self, event):
        """键盘按键事件的处理函数，用于实现撤回功能。"""
        if event.key == 'z':
            # 如果当前正在定义一条线（已点击了一个点）
            if self.markers:
                # 移除记录的点和对应的标记
                self.points.pop()
                last_marker = self.markers.pop()
                last_marker.remove()
                self.ax.figure.canvas.draw()
                print("上一个点已撤回。")
            # 如果要撤回一条已完成的线
            elif self.lines:
                # 移除上一条绘制的直线
                last_line = self.lines.pop()
                last_line.remove()
                self.ax.legend() # 更新图例
                self.ax.figure.canvas.draw()
                print("上一条线已撤回。")
            else:
                print("没有可撤回的操作。")

    def draw_line(self):
        """根据记录的两个点计算并绘制直线。"""
        x1, y1 = self.points[0]
        x2, y2 = self.points[1]

        # 计算直线方程 y = kx + b
        if abs(x2 - x1) < 1e-6: # 处理斜率不存在的垂直线情况
            k, b, eq_text = np.inf, np.nan, f"x = {x1:.2f}"
        else:
            k = (y2 - y1) / (x2 - x1)
            b = y1 - k * x1
            eq_text = f"y = {k:.4f}x + {b:.4f}"

        # 根据坐标轴的当前范围，计算直线的端点并绘制
        x_vals = np.array(self.ax.get_xlim())
        if np.isinf(k): # 绘制垂直线
             y_vals = np.array(self.ax.get_ylim())
             line = self.ax.plot([x1, x1], y_vals, 'purple', linestyle='--', linewidth=2, label=f'Interactive: {eq_text}')
        else: # 绘制普通斜率的直线
            y_vals = k * x_vals + b
            line = self.ax.plot(x_vals, y_vals, 'purple', linestyle='--', linewidth=2, label=f'Interactive: {eq_text}')
        
        self.lines.extend(line)
        
        # 清理本次操作的临时标记点和坐标
        for marker in self.markers:
            marker.remove()
        self.markers.clear()
        self.points.clear()

        # 更新图例并重绘图像
        self.ax.legend()
        self.ax.figure.canvas.draw()
        print("\n准备就绪，可继续点击绘制下一条直线。")

def plot_track_area_horizontal_glorys(DS: list, no: int, needed_idx: int, variable: str = 'vorticity',
                                   show_fig: bool = False, save_fig: bool = False, deep_argo: bool = False,
                                   k: float | list[float] | None = None, b: float | list[float] | None = None, 
                                   needed_depth: float | int = 0, inline_mode: bool = True):
    '''
    绘制指定涡旋在特定时刻的表层物理场快照及相关的Argo浮标数据。

    该函数支持两种模式：
    1. inline_mode=True (默认): 适用于静态图表生成和高分辨率保存，行为与原始版本完全一致。
    2. inline_mode=False: 适用于在Jupyter Notebook中使用 %matplotlib widget 进行交互式分析。

    参数:
        DS (list): 包含所有涡旋轨迹信息的数据集。
        no (int): 需要绘制的涡旋的唯一编号。
        needed_idx (int): 涡旋轨迹的时间点索引，用于确定绘图的具体日期。
        variable (str): 作为背景场绘制的GLORYS物理变量。默认为 'vorticity'。
        show_fig (bool): 是否在运行时显示生成的图像。默认为 False。
        save_fig (bool): 是否将生成的图像保存为文件。默认为 False。
        deep_argo (bool): 是否使用深层Argo数据模式。默认为 False。
        k (float | list[float], optional): 直线方程 y=kx+b 的斜率或斜率列表。默认为 None。
        b (float | list[float], optional): 直线方程 y=kx+b 的截距或截距列表。默认为 None。
        needed_depth (float | int): 需要绘制的GLORYS数据深度，默认为0（表层）。
        inline_mode (bool): 是否为静态内联模式。默认为True。设为False以启用交互式widget模式的优化。
    '''
    # 根据模式定义一套协调的尺寸参数
    if inline_mode:
        # 适用于高分辨率保存的静态模式尺寸，与原始版本完全一致
        figsize = (30, 25)
        title_fs, label_fs, tick_fs, legend_fs = 20, 20, 16, 18
        cbar_label_fs, cbar_tick_fs = 20, 14
        argo_text_fs = 7
        track_lw, contour_lw, circle_lw, line_lw = 1.0, 1.0, 1.0, 2.0
        cbar_pad = 0.046 # 原始cbar间距
    else:
        # 适用于交互式widget的屏幕友好尺寸
        figsize = (12, 10) 
        title_fs, label_fs, tick_fs, legend_fs = 16, 14, 12, 10
        cbar_label_fs, cbar_tick_fs = 12, 10
        argo_text_fs = 6
        track_lw, contour_lw, circle_lw, line_lw = 1.0, 1.0, 1.0, 2.0
        cbar_pad = 0.12 # 交互模式下增大cbar间距

    wanted_track = find_track(DS, no)
    num, time, center_lon, center_lat, _, _, contour_lon, contour_lat, radius, _, _ = zip(*wanted_track)
    dates = convert_date(time) if time else None

    # 获取Argo浮标数据
    argo_data_filtered = filtered_float_data(DS, no)
    argo_data_filtered = argo_data_filtered[pd.to_datetime(argo_data_filtered[['Year', 'Month', 'Day']])==dates[needed_idx]]
    
    if deep_argo:
        filtered_by_depth = argo_data_filtered[argo_data_filtered['Depth'] >= 500].copy()
        if filtered_by_depth.empty:
            print("Warning: No data found with Depth >= 500.")
            needed_data = pd.DataFrame(columns=argo_data_filtered.columns)
        else:
            idx_max_do = filtered_by_depth.groupby('Profile_number')['DO'].idxmax()
            needed_data = filtered_by_depth.loc[idx_max_do]
            needed_data.index.name = None
    else:
        if argo_data_filtered.empty:
            needed_data = pd.DataFrame(columns=argo_data_filtered.columns)
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
        glorys_lon_filtered, glorys_lat_filtered, glorys_depth_filtered, glorys_variables_filtered = get_track_area_glorys(DS, no, needed_idx, variables=['u', 'v'], depth=needed_depth)
        zeta, f = calculate_vorticity(glorys_lon_filtered, glorys_lat_filtered, glorys_variables_filtered['u'], glorys_variables_filtered['v'])
        glorys_variable_filtered = zeta/f
    else:
        glorys_lon_filtered, glorys_lat_filtered, glorys_depth_filtered, glorys_variables_filtered = get_track_area_glorys(DS, no, needed_idx, variables=[variable], depth=needed_depth)
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

    fig, ax = plt.subplots(figsize=figsize)
    ax.set_title(f'Track {ds_names}{num[0]} at {glorys_depth_filtered[0]:.2f}m, {dates[needed_idx].strftime("%Y-%m-%d")}', fontsize=title_fs)
    ax.set_xlabel('Longitude', fontsize=label_fs)
    ax.set_ylabel('Latitude', fontsize=label_fs)
    world.plot(color='green', ax=ax)

    ax.tick_params(axis='both', which='major', labelsize=tick_fs)

    # 绘制涡旋轨迹
    ax.plot(center_lon, center_lat, color=colors, linewidth=track_lw, label='Center Track')
    ax.plot(center_lon[0], center_lat[0], marker='o', color=colors, markersize=10)
    ax.plot(center_lon[-1], center_lat[-1], marker='x', color=colors, markersize=10)

    # 绘制背景场
    pc = ax.pcolormesh(glorys_lon_filtered, glorys_lat_filtered, glorys_variable_filtered, cmap='seismic', shading='auto', alpha=0.5)
    cbar = plt.colorbar(pc, ax=ax, orientation='horizontal', fraction=0.046, pad=0.046)
    if variable == 'vorticity':
        cbar.set_label(r'$\zeta/f$', fontsize=cbar_label_fs)
        pc.set_clim(-0.7,0.7)
    elif variable == 'thetao':
        cbar.set_label('Temperature (°C)', fontsize=cbar_label_fs)
    elif variable == 'so':
        cbar.set_label('Salinity (psu)', fontsize=cbar_label_fs)
    elif variable == 'u':
        cbar.set_label('Zonal Velocity (m/s)', fontsize=cbar_label_fs)
    elif variable == 'v':
        cbar.set_label('Meridional Velocity (m/s)', fontsize=cbar_label_fs)
    elif variable == 'ssh':
        cbar.set_label('Sea Surface Height (m)', fontsize=cbar_label_fs)
    else:
        cbar.set_label(variable, fontsize=cbar_label_fs)
    cbar.ax.tick_params(labelsize=cbar_tick_fs)

    # 绘制Argo浮标数据
    if not needed_data.empty:
        if deep_argo:
            sc = ax.scatter(needed_data['Longitude'], needed_data['Latitude'], c=needed_data['DO'], cmap = 'bwr', s=180,
                            vmin=150, vmax=240, edgecolors='black', linewidths=0.5, label='Argo with max DO under 500m', zorder=5)
            cbar2 = plt.colorbar(sc, ax=ax, orientation='horizontal', fraction=0.046, pad=cbar_pad)
            cbar2.set_label('DO/μmol·kg⁻¹', fontsize=cbar_label_fs)
            cbar2.ax.tick_params(labelsize=cbar_tick_fs)
        else:
            ax.scatter(needed_data['Longitude'], needed_data['Latitude'], color='blue', s=180, label='Argo', zorder=5)
        for idx, row in needed_data.iterrows():
            ax.text(row['Longitude'], row['Latitude'], f"{int(row['Depth'])}", fontsize=argo_text_fs, fontweight='bold', ha='center', va='center', color='black', zorder=6)
    else:
        print(f"No Argo data available for eddy {ds_names}{no} at the specified index {needed_idx}.")

    # 绘制当前时刻涡旋
    circle = plt.Circle((center_lon[needed_idx], center_lat[needed_idx]), radius[needed_idx] / 111320.0,
                        color='r', fill=False, linestyle='--', alpha=0.2, linewidth=circle_lw, label='Effective Radius')
    ax.add_patch(circle)
    ax.scatter(center_lon[needed_idx], center_lat[needed_idx], color='black', s=20, label='Eddy Center', zorder=5)
    ax.plot(contour_lon[needed_idx], contour_lat[needed_idx], color=colors, linewidth=contour_lw, alpha=0.5, label='Effective Contour')

    # 绘制 y = kx + b 直线
    if k is not None and b is not None:
        k_list = [k] if isinstance(k, (int, float)) else k
        b_list = [b] if isinstance(b, (int, float)) else b

        if len(k_list) != len(b_list):
            raise ValueError("The lists for k and b must have the same length.")

        line_x = np.array([glorys_lon_min, glorys_lon_max])

        for i, (k_val, b_val) in enumerate(zip(k_list, b_list)):
            line_y = k_val * line_x + b_val
            ax.plot(line_x, line_y, color='purple', linestyle='-', linewidth=line_lw, label=f'Profile Line {i+1}: y={k_val:.2f}x+{b_val:.2f}')

    ax.legend(fontsize=legend_fs)
    ax.set_xlim(glorys_lon_min, glorys_lon_max)
    ax.set_ylim(glorys_lat_min, glorys_lat_max)
    ax.set_aspect('equal')

    # 紧凑布局，消除多余空白
    plt.tight_layout()

    # 保存图片
    if save_fig:
        output_dir = "plot_track_area_horizontal_glorys"
        os.makedirs(output_dir, exist_ok=True)
        base_filename = f"{ds_names}{no}_{glorys_depth_filtered[0]:.2f}m_{variable}_{dates[needed_idx].strftime('%Y%m%d')}.png"
        plt.savefig(os.path.join(output_dir, base_filename), dpi=300, bbox_inches='tight')
        print(f"\nFigure saved to: {os.path.join(output_dir, base_filename)}")

    # 显示图像
    if show_fig:
        # 仅在非内联模式（即交互模式）下，激活直线绘制器
        if not inline_mode:
            line_drawer = LineDrawer(ax)
            fig.canvas.mpl_connect('button_press_event', line_drawer.onclick)
        
        plt.show()

    # 只有在静态内联模式下，才在函数结束时关闭图像以释放内存
    if inline_mode:
        plt.close(fig)

def get_track_area_glorys(DS: list, no: int, needed_idx: int | pd.Timestamp, variables: list = ['thetao'], depth: float | int | None = None):
    '''
    获取指定涡旋在特定时间点周围的 GLORYS 数据。

    该函数会根据涡旋轮廓确定一个矩形区域，并从相应的 GLORYS 文件中
    提取此区域内的一个或多个物理变量。

    参数:
        DS (list): 包含涡旋轨迹信息的数据集。
        no (int): 涡旋的唯一编号。
        needed_idx (int | pd.Timestamp): 需要提取数据的时间点索引或时间戳。
        variables (list): 需要提取的变量列表，默认为 ['thetao']，可选'salinity', 'u', 'v', 'ssh', 'mlt'。
        depth (float | int | None): 如果指定，提取该深度的 GLORYS 数据；如果为 None，则提取2000米以内的所有深度数据。

    返回:
        一个元组，包含筛选后的经度、纬度、深度数组，以及一个存储了所有请求变量数据的字典。
    '''
    wanted_track = find_track(DS, no)
    _, _, _, _, _, _, contour_lon, contour_lat, _, _, _ = zip(*wanted_track)

    if type(needed_idx) is int:
        glorys_filepaths_dict = find_track_glorys_filepath(DS, no)
        needed_glorys_data = Dataset(list(glorys_filepaths_dict.values())[needed_idx], 'r')
    elif isinstance(needed_idx, pd.Timestamp):
        glorys_filepaths_dict = {needed_idx: get_glorys_filepath(needed_idx)}
        needed_glorys_data = Dataset(glorys_filepaths_dict[needed_idx], 'r')
    else:
        raise ValueError("needed_idx must be an integer index or a pd.Timestamp.")

    contour_lon_filtered = np.ma.masked_equal(contour_lon, 180.0)
    contour_lat_filtered = np.ma.masked_equal(contour_lat, 0.0)

    glorys_lon_min = np.min(contour_lon_filtered) - 0.5
    glorys_lon_max = np.max(contour_lon_filtered) + 0.5
    glorys_lat_min = np.min(contour_lat_filtered) - 0.5
    glorys_lat_max = np.max(contour_lat_filtered) + 0.5

    glorys_lon = needed_glorys_data.variables['longitude'][:]
    glorys_lat = needed_glorys_data.variables['latitude'][:]
    glorys_depth = needed_glorys_data.variables['depth'][:]

    glorys_lon_mask = (glorys_lon >= glorys_lon_min) & (glorys_lon <= glorys_lon_max)
    glorys_lat_mask = (glorys_lat >= glorys_lat_min) & (glorys_lat <= glorys_lat_max)
    if depth is not None:
        glorys_depth_mask = np.zeros_like(glorys_depth, dtype=bool)
        if glorys_depth.size > 0:
            k = np.argmin(np.abs(glorys_depth - depth))
            glorys_depth_mask[k] = True
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
            if depth is not None:
                glorys_variables_filtered['thetao'] = glorys_variables_filtered['thetao'][0, :, :]
        elif var == 'salinity' or var == 'so':
            glorys_variables_filtered['salinity'] = needed_glorys_data.variables['so'][0, glorys_depth_mask, glorys_lat_mask, glorys_lon_mask]
            if depth is not None:
                glorys_variables_filtered['salinity'] = glorys_variables_filtered['salinity'][0, :, :]
        elif var == 'u' or var == 'uo':
            glorys_variables_filtered['u'] = needed_glorys_data.variables['uo'][0, glorys_depth_mask, glorys_lat_mask, glorys_lon_mask]
            if depth is not None:
                glorys_variables_filtered['u'] = glorys_variables_filtered['u'][0, :, :]
        elif var == 'v' or var == 'vo':
            glorys_variables_filtered['v'] = needed_glorys_data.variables['vo'][0, glorys_depth_mask, glorys_lat_mask, glorys_lon_mask]
            if depth is not None:
                glorys_variables_filtered['v'] = glorys_variables_filtered['v'][0, :, :]
        elif var == 'ssh' or var == 'zos':
            glorys_variables_filtered['ssh'] = needed_glorys_data.variables['zos'][0, glorys_lat_mask, glorys_lon_mask]
        elif var == 'mlt' or var == 'mlotst':
            glorys_variables_filtered['mlt'] = needed_glorys_data.variables['mlotst'][0, glorys_lat_mask, glorys_lon_mask]
        else:
            raise ValueError(f"Unsupported variable: {var}. Supported variables are: 'thetao', 'so', 'uo', 'vo', 'zos', 'mlotst'.")

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
    
def get_idx(DS: list, no: int, start_date: str, end_date: str = None) -> int | list | None:
    '''
    获取指定涡旋编号在给定时间或时间范围内的索引。

    - 如果只提供 start_date，则返回该日期对应的单个索引。
    - 如果同时提供 start_date 和 end_date，则返回该时间范围内的索引列表。

    参数:
        DS (list): 涡旋轨迹数据集。
        no (int): 涡旋编号。
        start_date (str): 起始日期，格式为 'YYYY-MM-DD'。
        end_date (str, optional): 结束日期，格式为 'YYYY-MM-DD'。默认为 None。

    返回:
        int | list | None: 
        - 如果 end_date 为 None，返回单个整数索引或 None (如果未找到)。
        - 如果提供了 end_date，返回一个整数索引的列表。
    '''
    wanted_track = find_track(DS, no)
    if not wanted_track:
        print(f"未找到涡旋 {no} 的轨迹数据。")
        # 根据调用模式返回正确的空值
        return None if end_date is None else []

    start_date_dt = pd.to_datetime(start_date)

    # 情况一：只查找单个日期的索引
    if end_date is None:
        for idx, track_point in enumerate(wanted_track):
            track_date = convert_date(track_point[1])
            # 精确匹配日期，忽略时间部分
            if track_date.date() == start_date_dt.date():
                return idx  # 找到后立即返回单个索引
        
        # 如果循环结束仍未找到匹配的日期
        print(f"在涡旋 {no} 的轨迹中未找到日期 {start_date}。")
        return None

    # 情况二：查找一个日期范围内的索引列表
    else:
        end_date_dt = pd.to_datetime(end_date)
        idx_list = []
        for idx, track_point in enumerate(wanted_track):
            track_date = convert_date(track_point[1])
            if start_date_dt <= track_date <= end_date_dt:
                idx_list.append(idx)
        
        return idx_list

def get_vertical_glorys(DS: list, no: int, needed_idx: int,
                        k: float | list[float], b: float | list[float],
                        variables: list = ['vorticity']) -> list[dict]:
    '''
    计算并返回指定涡旋在特定时刻沿一条或多条 y=kx+b 剖面的物理量数据。

    该函数封装了从GLORYS数据场中提取剖面数据的核心插值计算，并以字典列表的形式返回结果。
    它能正确处理三维变量（如温度，返回二维垂直剖面）和二维变量（如混合层深度，返回一维水平剖面）。
    无论输入变量使用何种别名，输出字典中的键都将是标准化的变量名。

    参数:
        DS (list): 包含涡旋轨迹信息的数据集。
        no (int): 涡旋的唯一编号。
        needed_idx (int): 涡旋轨迹的时间点索引。
        k (float | list[float]): 直线方程 y=kx+b 的斜率或斜率列表。
        b (float | list[float]): 直线方程 y=kx+b 的截距或截距列表。
        variables (list): 需要提取的GLORYS物理变量列表。默认为 ['vorticity']。

    返回:
        list[dict]: 一个包含一个或多个剖面结果字典的列表。列表中每个字典的结构如下：
        
          - **'profile_data' (dict)**:
            - *键*: 标准化的物理量名称 (字符串, 如 'salinity', 'mlt')。
            - *值*:
                - 对于三维变量 (如'vorticity', 'thetao'), 值为一个二维 `numpy.ma.MaskedArray` 数组，代表其垂直剖面，维度为 `(深度层数, 剖面水平点数)`。
                - 对于二维变量 (如'mlt', 'ssh'), 值为一个一维 `numpy.ndarray` 数组，代表其沿剖面线的水平分布，长度为 `剖面水平点数`。

          - **'y_coords' (np.ndarray)**:
            - 一个一维NumPy数组，表示剖面的横坐标轴（物理距离）。
            - 数值单位为公里 (km)，`0` 点对应涡旋中心在剖面线上的投影位置。

          - **'z_coords' (np.ndarray)**:
            - 一个一维NumPy数组，表示剖面的纵坐标轴（物理深度）。
            - 数值单位为米 (m)，代表了GLORYS数据中的深度分层。
          
          - **'lon_coords' (np.ndarray)**:
            - 一个一维NumPy数组，包含了剖面线上每个点（对应`y_coords`）的经度。
          
          - **'lat_coords' (np.ndarray)**:
            - 一个一维NumPy数组，包含了剖面线上每个点（对应`y_coords`）的纬度。

          - **'projections' (dict)**:
            - 一个字典，包含了涡旋边界在横坐标 (`y_coords`) 上的投影位置。
            - *键*: 边界类型 (字符串)， `'radius'` 代表有效半径，`'contour'` 代表有效轮廓。
            - *值*: 一个包含交点位置(km)的列表。

          - **'metadata' (dict)**:
            - 一个包含涡旋元数据的字典。
            - *键*:
                - 'eddy_no'` (int): 涡旋的唯一编号。
                - 'date_str'` (str): 该剖面对应日期的字符串，格式为 'YYYY-MM-DD'。
                - 'k'` (float): 该剖面所用直线方程的斜率。
                - 'b'` (float): 该剖面所用直线方程的截距。
    '''
    # --- 0. 准备工作：统一输入格式并获取公共数据 ---
    if k is None or b is None:
        raise ValueError("k 和 b 必须提供以计算垂直剖面。")

    k_list = [k] if isinstance(k, (int, float)) else k
    b_list = [b] if isinstance(b, (int, float)) else b

    if len(k_list) != len(b_list):
        raise ValueError("k 和 b 的列表长度必须一致。")

    # 建立别名到标准名，以及变量维度的映射
    alias_map = {
        'thetao': 'thetao',
        'salinity': 'salinity', 'so': 'salinity',
        'u': 'u', 'uo': 'u',
        'v': 'v', 'vo': 'v',
        'ssh': 'ssh', 'zos': 'ssh',
        'mlt': 'mlt', 'mlotst': 'mlt',
        'vorticity': 'vorticity'
    }
    var_dims = {
        'thetao': 3, 'salinity': 3, 'u': 3, 'v': 3, 'vorticity': 3,
        'ssh': 2, 'mlt': 2
    }

    raw_vars_to_fetch = set()
    for var in variables:
        if var == 'vorticity': raw_vars_to_fetch.update(['u', 'v'])
        else: raw_vars_to_fetch.add(var)

    wanted_track = find_track(DS, no)
    if not wanted_track:
        print(f"未找到涡旋 {no} 的轨迹数据。")
        return [{} for _ in k_list]

    num, time, center_lon, center_lat, _, _, contour_lon, contour_lat, radius, _, _ = zip(*wanted_track)
    dates = convert_date(time)
    
    glorys_lon_raw, glorys_lat_raw, glorys_depth_raw, glorys_data_raw = get_track_area_glorys(
        DS, no, needed_idx, variables=list(raw_vars_to_fetch)
    )
    
    if glorys_depth_raw.size == 0 and not all(var_dims.get(alias_map.get(v, v)) == 2 for v in variables):
        return [{} for _ in k_list]

    all_profiles_data = []

    # --- 开始循环，为每一对 k, b 计算一个剖面 ---
    for k_val, b_val in zip(k_list, b_list):
        # --- 1. 计算水平剖面线的坐标 ---
        R_earth = 6371e3
        contour_lon_filtered = np.ma.masked_equal(contour_lon, 180.0)
        glorys_lon_min, glorys_lon_max = np.min(contour_lon_filtered) - 0.5, np.max(contour_lon_filtered) + 0.5
        contour_lat_filtered = np.ma.masked_equal(contour_lat, 0.0)
        glorys_lat_min, glorys_lat_max = np.min(contour_lat_filtered) - 0.5, np.max(contour_lat_filtered) + 0.5
        
        num_points = 500
        if k_val == 0:
            profile_lons = np.linspace(glorys_lon_min, glorys_lon_max, num_points)
            profile_lats = np.full_like(profile_lons, b_val)
        else:
            lons_temp = np.linspace(glorys_lon_min, glorys_lon_max, num_points)
            lats_temp = k_val * lons_temp + b_val
            mask = (lats_temp >= glorys_lat_min) & (lats_temp <= glorys_lat_max)
            profile_lons, profile_lats = lons_temp[mask], lats_temp[mask]
            if len(profile_lons) < 2: 
                all_profiles_data.append({})
                continue

        dlat = np.deg2rad(np.diff(profile_lats)); dlon = np.deg2rad(np.diff(profile_lons))
        mid_lats = np.deg2rad((profile_lats[:-1] + profile_lats[1:]) / 2)
        dist_segments = R_earth * np.sqrt(dlat**2 + (np.cos(mid_lats) * dlon)**2)
        y_coords_raw = np.insert(np.cumsum(dist_segments), 0, 0) / 1000.0
        
        current_center_lon, current_center_lat = center_lon[needed_idx], center_lat[needed_idx]
        if k_val == 0: xp, yp = current_center_lon, b_val
        else:
            xp = (current_center_lon + k_val * current_center_lat - k_val * b_val) / (1 + k_val**2)
            yp = k_val * xp + b_val
        center_idx_on_profile = np.argmin((profile_lons - xp)**2 + (profile_lats - yp)**2)
        y_coords_recenter = y_coords_raw - y_coords_raw[center_idx_on_profile]

        # --- 2. 区分维度，计算剖面数据 ---
        z_coords = glorys_depth_raw
        profile_data_dict = {}
        for var in variables:
            standard_name = alias_map.get(var, var)
            dimension = var_dims.get(standard_name, 3)

            # --- 2.1 处理2D变量 (如 mlt, ssh) ---
            if dimension == 2:
                data_2d = glorys_data_raw.get(standard_name)
                if data_2d is None or data_2d.size == 0 or np.all(np.ma.getmask(data_2d)):
                    profile_data_dict[standard_name] = np.ma.masked_all(len(profile_lons))
                    continue
                
                interp_func_2d = RegularGridInterpolator((glorys_lat_raw, glorys_lon_raw), data_2d.filled(np.nan), 
                                                         method='linear', bounds_error=False, fill_value=np.nan)
                interpolated_1d = interp_func_2d(list(zip(profile_lats, profile_lons)))
                profile_data_dict[standard_name] = np.ma.masked_invalid(interpolated_1d)
                continue

            # --- 2.2 处理3D变量 (如 thetao, vorticity) ---
            glorys_variable_3d = None
            if standard_name == 'vorticity':
                u, v = glorys_data_raw.get('u'), glorys_data_raw.get('v')
                if u is not None and v is not None and u.size > 0 and v.size > 0:
                    if u.ndim == 2: u, v = u[np.newaxis, :, :], v[np.newaxis, :, :]
                    zeta_3d, f_3d = calculate_vorticity(glorys_lon_raw, glorys_lat_raw, u, v)
                    glorys_variable_3d = zeta_3d / f_3d
            else:
                glorys_variable_3d = glorys_data_raw.get(standard_name)

            if glorys_variable_3d is None or glorys_variable_3d.size == 0 or np.all(np.ma.getmask(glorys_variable_3d)):
                profile_data_dict[standard_name] = np.ma.masked_all((len(z_coords), len(profile_lons)))
                continue

            if glorys_variable_3d.ndim == 2: glorys_variable_3d = glorys_variable_3d[np.newaxis, :, :]
            
            query_depths, query_lats = np.meshgrid(z_coords, profile_lats, indexing='ij')
            _, query_lons = np.meshgrid(z_coords, profile_lons, indexing='ij')
            xi_points = np.vstack([query_depths.ravel(), query_lats.ravel(), query_lons.ravel()]).T
            
            interp_func = RegularGridInterpolator((z_coords, glorys_lat_raw, glorys_lon_raw), glorys_variable_3d.filled(np.nan), method='linear', bounds_error=False, fill_value=np.nan)
            interpolated_values_flat = interp_func(xi_points)
            profile_data_dict[standard_name] = np.ma.masked_invalid(interpolated_values_flat.reshape(len(z_coords), len(profile_lons)))

        # --- 3. 计算边界投影 ---
        effective_radius_deg = radius[needed_idx] / 111320.0
        A, B = 1 + k_val**2, 2 * (k_val*b_val - k_val*current_center_lat - current_center_lon)
        C = current_center_lon**2 + (b_val - current_center_lat)**2 - effective_radius_deg**2
        discriminant = B**2 - 4*A*C
        radius_intersections_lon = [(-B + s * np.sqrt(discriminant)) / (2*A) for s in [-1, 1]] if discriminant >= 0 else []
        radius_proj_dists = [y_coords_raw[np.argmin((profile_lons - lon_i)**2 + (profile_lats - (k_val*lon_i + b_val))**2)] - y_coords_raw[center_idx_on_profile] for lon_i in radius_intersections_lon]

        contour_lon_valid = contour_lon[needed_idx][contour_lon[needed_idx] != 180.0]
        contour_lat_valid = contour_lat[needed_idx][contour_lat[needed_idx] != 0.0]
        contour_intersections_xy = find_polygon_line_intersections(contour_lon_valid, contour_lat_valid, profile_lons, profile_lats)
        contour_proj_dists = [y_coords_raw[np.argmin((profile_lons - lon_i)**2 + (profile_lats - lat_i)**2)] - y_coords_raw[center_idx_on_profile] for lon_i, lat_i in contour_intersections_xy]

        projections_dict = {'radius': sorted(radius_proj_dists), 'contour': sorted(contour_proj_dists)}
        
        # --- 4. 整合单个剖面的结果 ---
        single_profile_result = {
            'profile_data': profile_data_dict,
            'y_coords': y_coords_recenter,
            'z_coords': z_coords,
            'lon_coords': profile_lons,
            'lat_coords': profile_lats,
            'projections': projections_dict,
            'metadata': {
                'eddy_no': num[0],
                'date_str': dates[needed_idx].strftime("%Y-%m-%d"),
                'k': k_val,
                'b': b_val,
            }
        }
        all_profiles_data.append(single_profile_result)

    # --- 5. 返回所有剖面的结果列表 ---
    return all_profiles_data

def plot_vertical_glorys(DS: list, no: int, needed_idx: int, 
                         k: float | list[float], b: float | list[float], 
                         variable: str = 'vorticity',
                         show_fig: bool = False, save_fig: bool = False, 
                         xmin: float = None, xmax: float = None,
                         ymin: float = None, ymax: float = None,
                         plot_mlt: bool = False):
    '''
    绘制指定涡旋在特定时刻，沿一条或多条剖面线 (y = kx + b) 的物理量垂直剖面图。

    该函数调用 get_vertical_glorys 获取指定变量的数据，然后进行可视化。
    当 k 和 b 为列表时，会为每一对 (k, b) 生成一张独立的图表。
    可以选择性地在图上叠加混合层深度曲线。

    参数:
        (同 get_vertical_glorys, 但 variable 为 str，一次只处理一个变量)
        ...
        show_fig (bool): 是否显示图像。
        save_fig (bool): 是否保存图像。
        xmin, xmax, ymin, ymax (float): 坐标轴范围。
        plot_mlt (bool): 是否在图上绘制混合层深度分界线。默认为 False。
    '''
    # --- 1. 获取所有计算好的剖面数据包 ---
    
    vars_to_fetch = {variable}
    if plot_mlt:
        vars_to_fetch.add('mlt')
        
    all_data_packages = get_vertical_glorys(DS, no, needed_idx, k, b, variables=list(vars_to_fetch))

    if not all_data_packages:
        print(f"警告: get_vertical_glorys 未能返回任何数据。绘图已取消。")
        return

    # --- 开始循环，为每一个数据包生成一张图 ---
    for data_package in all_data_packages:
        if not data_package:
            print(f"警告: 收到一个空的数据包，跳过此剖面的绘图。")
            continue
            
        profile_variable_2d = data_package['profile_data'].get(variable)
        if profile_variable_2d is None:
             alias_map = {'so': 'salinity', 'uo': 'u', 'vo': 'v'}
             standard_name = alias_map.get(variable, variable)
             profile_variable_2d = data_package['profile_data'].get(standard_name)

        if profile_variable_2d is None or np.all(getattr(profile_variable_2d, 'mask', True)):
            k_meta, b_meta = data_package.get('metadata', {}).get('k'), data_package.get('metadata', {}).get('b')
            print(f"警告: 变量 '{variable}' 在剖面 k={k_meta}, b={b_meta} 上的数据无效。绘图已取消。")
            continue

        # --- 2. 准备绘图 ---
        y_coords = data_package['y_coords']
        z_coords = data_package['z_coords']
        projections = data_package['projections']
        metadata = data_package['metadata'] 

        callers_local_vars = inspect.currentframe().f_back.f_locals
        ds_name = [var_name for var_name, var_val in callers_local_vars.items() if var_val is DS][0].upper()
        
        prop_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        eddy_color = prop_colors[1] if 'AC' in ds_name else prop_colors[0]
        
        if profile_variable_2d.shape[1] != len(y_coords):
           y_coords = y_coords[:profile_variable_2d.shape[1]]

        Y_mesh, Z_mesh = np.meshgrid(y_coords, z_coords)
        
        # 设置变量相关的绘图属性
        if variable == 'vorticity': cbar_label, cmap, clim = r'$\zeta/f$', 'seismic', (-0.3, 0.3)
        elif variable in ['thetao']: cbar_label, cmap = 'Temperature (°C)', 'rainbow'
        elif variable in ['salinity', 'so']: cbar_label, cmap = 'Salinity (psu)', 'viridis'
        elif variable in ['u', 'v', 'uo', 'vo']: cbar_label, cmap = 'Velocity (m/s)', 'RdBu_r'
        else: cbar_label, cmap, clim = variable, 'viridis', None
        
        if 'clim' not in locals() or clim is None:
            valid_values = profile_variable_2d[~profile_variable_2d.mask]
            clim = (valid_values.min(), valid_values.max()) if valid_values.size > 0 else (0,1)
            if variable in ['u', 'v', 'uo', 'vo']:
                max_abs = np.max(np.abs(valid_values)) if valid_values.size > 0 else 1
                clim = (-max_abs, max_abs)

        # --- 3. 执行绘图 ---
        fig, ax = plt.subplots(figsize=(20, 15))
        
        date_str = metadata['date_str']
        title = (f"Vertical Profile of {cbar_label} for Track {ds_name}{metadata['eddy_no']} "
                 f"on {date_str}, y={metadata['k']:.2f}x+{metadata['b']:.2f}")
        ax.set_title(title, fontsize=20)
        ax.set_xlabel('Distance from Eddy Center Projection (km)', fontsize=18)
        ax.set_ylabel('Depth (m)', fontsize=18)
        ax.tick_params(labelsize=14)
        
        pc = ax.pcolormesh(Y_mesh, Z_mesh, profile_variable_2d, cmap=cmap, shading='auto', vmin=clim[0], vmax=clim[1])
        cbar = fig.colorbar(pc, ax=ax, orientation='vertical', fraction=0.046, pad=0.04)
        cbar.set_label(cbar_label, fontsize=18)
        cbar.ax.tick_params(labelsize=14)

        ax.axvline(0, color='black', linestyle='--', linewidth=2, label='Eddy Center Projection')
        for i, dist in enumerate(projections['radius']): 
            ax.axvline(dist, color='r', linestyle='--', linewidth=2, label='Effective Radius Projection' if i == 0 else "")
        for i, dist in enumerate(projections['contour']): 
            ax.axvline(dist, color=eddy_color, linestyle=':', linewidth=2, label='Effective Contour Projection' if i == 0 else "")
            
        # 绘制混合层深度
        mld_lines = () # 创建一个空元组，用于存放MLD的两条线
        if plot_mlt:
            mlt_data = data_package['profile_data'].get('mlt')
            if mlt_data is not None and not np.all(getattr(mlt_data, 'mask', True)):
                # **捕获两条线的艺术家对象**
                # 注意 plot 返回的是一个列表，所以我们用 l, 来解包
                black_line, = ax.plot(y_coords, mlt_data, color='black', linewidth=2.5, zorder=3)
                white_line, = ax.plot(y_coords, mlt_data, color='white', linewidth=1.5, zorder=4, label='Mixed Layer Depth')
                mld_lines = (black_line, white_line) # 将两条线存入元组
            else:
                print(f"注意: 未能在剖面 k={metadata['k']}, b={metadata['b']} 上找到有效的混合层深度数据。")

        ax.set_ylim(z_coords.max(), z_coords.min())
        if xmin is not None and xmax is not None: ax.set_xlim(xmin, xmax)
        if ymin is not None and ymax is not None: ax.set_ylim(ymax, ymin)
        
        # --- 构建并自定义图例 ---
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles)) # 创建一个去重的标签-句柄字典

        # 如果我们画了MLD线，就用我们创建的复合句柄替换掉自动生成的单个白线句柄
        if plot_mlt and mld_lines:
            by_label['Mixed Layer Depth'] = mld_lines
            
        ax.legend(by_label.values(), by_label.keys(), fontsize=18)

        # --- 4. 保存和显示 ---
        if save_fig:
            output_dir = "plot_vertical_glorys"
            os.makedirs(output_dir, exist_ok=True)
            date_fn = date_str.replace('-', '')
            base_filename = (f"{ds_name}{metadata['eddy_no']}_vertical_{variable}_{date_fn}_"
                             f"k{metadata['k']:.2f}b{metadata['b']:.2f}.png")
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
    if o4 == 0 and _on_segment(p2, q1, q2): return q1

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
    
    if len(polygon_lon) < 3: # 需要至少三个点才能构成多边形
        return []

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

def regrid_vertical_slice(data_package: dict, dy: float = None, dz: float = None) -> dict:
    '''
    将垂直剖面数据包插值到新的等间距网格上。

    该函数接收 get_vertical_glorys 的输出，并根据用户指定的水平(dy)或
    垂直(dz)间距，生成一个新的、网格化更规整的数据包。

    参数:
        data_package (dict): 从 get_vertical_glorys 函数获取的原始数据包。
        dy (float, optional): 新的水平(y轴)网格间距，单位为公里(km)。
                              默认为 None，表示不改变水平网格。
        dz (float, optional): 新的垂直(z轴)网格间距，单位为米(m)。
                              默认为 None，表示不改变垂直网格。

    返回:
        dict: 一个结构与输入相同但数据已被插值到新网格上的新数据包。
              字典包含以下键值对：

            - **'profile_data' (dict)**:
                - 其内部结构与输入相同，但每个二维数组的值都是插值后的结果。

            - **'y_coords' (np.ndarray)**:
                - 如果提供了 `dy`，这将是一个新生成的一维等间距数组。
                - 否则，与输入的 'y_coords' 相同。

            - **'z_coords' (np.ndarray)**:
                - 如果提供了 `dz`，这将是一个新生成的一维等间距数组。
                - 否则，与输入的 'z_coords' 相同。
            
            - **'lon_coords' (np.ndarray)**:
                - 如果水平坐标被重采样，这将是新的一维插值经度数组。
                - 否则，与输入的 'lon_coords' 相同。

            - **'lat_coords' (np.ndarray)**:
                - 如果水平坐标被重采样，这将是新的一维插值纬度数组。
                - 否则，与输入的 'lat_coords' 相同。

            - **'projections' (dict)**:
                - 从原始数据包中原样复制而来，其数值仍对应原始的坐标系。
                
            - **'metadata' (dict)**:
                - 从原始数据包中原样复制而来。
    '''
    # 如果用户未指定任何新的网格间距，则直接返回原始数据包的深拷贝
    if dy is None and dz is None:
        return copy.deepcopy(data_package)

    # --- 1. 提取原始数据和坐标 ---
    original_y = data_package['y_coords']
    original_z = data_package['z_coords']
    original_profiles = data_package['profile_data']
    original_lon = data_package['lon_coords']
    original_lat = data_package['lat_coords']
    
    # RegularGridInterpolator 要求坐标必须是严格递增的。
    if np.any(np.diff(original_z) <= 0):
        raise ValueError("原始z坐标（深度）必须是严格递增的才能进行插值。")
    if np.any(np.diff(original_y) <= 0):
        raise ValueError("原始y坐标（距离）必须是严格递增的才能进行插值。")

    # --- 2. 创建新的网格坐标 ---
    new_y = np.arange(original_y.min(), original_y.max(), dy) if dy is not None else original_y
    if dz is not None:
        z_min, z_max = original_z.min(), original_z.max()
        new_z = np.arange(z_min, z_max, dz)
    else:
        new_z = original_z

    # --- 3. 对每个变量的剖面数据及经纬度进行插值 ---
    # 插值经纬度 (1D)，仅当水平坐标被重采样时执行
    if dy is not None:
        new_lon = np.interp(new_y, original_y, original_lon)
        new_lat = np.interp(new_y, original_y, original_lat)
    else:
        new_lon, new_lat = original_lon, original_lat

    # 插值剖面数据 (2D)
    new_profiles = {}
    # 创建查询点网格，用于高效插值, 'ij'索引确保维度顺序正确
    zz_grid, yy_grid = np.meshgrid(new_z, new_y, indexing='ij')
    query_points = np.vstack([zz_grid.ravel(), yy_grid.ravel()]).T

    for var, data_array in original_profiles.items():
        if data_array.mask.all():
            new_profiles[var] = np.ma.masked_all((len(new_z), len(new_y)))
            continue

        # 创建插值器实例
        interp_func = RegularGridInterpolator(
            (original_z, original_y), 
            data_array.filled(np.nan),
            method='linear', 
            bounds_error=False, 
            fill_value=np.nan
        )
        
        # 在新网格上执行插值
        new_data_flat = interp_func(query_points)
        
        # 将一维插值结果重塑为二维网格
        new_data = new_data_flat.reshape(len(new_z), len(new_y))
        
        # 将结果转回 masked array 并存入新字典
        new_profiles[var] = np.ma.masked_invalid(new_data)

    # --- 4. 构建并返回新的数据包 ---
    new_data_package = {
        'profile_data': new_profiles,
        'y_coords': new_y,
        'z_coords': new_z,
        'lon_coords': new_lon,
        'lat_coords': new_lat,
        'projections': copy.deepcopy(data_package['projections']),
        'metadata': copy.deepcopy(data_package['metadata'])
    }

    return new_data_package

def plot_data_package(data_package: dict, DS: list, variable: str,
                      show_fig: bool = False, save_fig: bool = False, xmin: float = None, xmax: float = None,
                      ymin: float = None, ymax: float = None):
    '''
    根据一个数据包和原始涡旋数据集，绘制单一物理量的垂直剖面图。

    该函数是一个灵活的可视化接口，它接收一个数据包，并根据传入的涡旋数据集
    (DS)来确定绘图风格（如颜色和标题），适用于 get_vertical_glorys 或
    regrid_vertical_slice 的输出。

    参数:
        data_package (dict): 从 get_vertical_glorys 或 regrid_vertical_slice 获取的数据包。
        DS (list): 原始的涡旋轨迹信息数据集 (如ACL, CL)，用于确定涡旋类型和名称。
        variable (str): 需要从 data_package 中绘制的变量名。
        show_fig (bool): 是否显示图像。
        save_fig (bool): 是否保存图像。
        xmin, xmax, ymin, ymax (float): 坐标轴范围。
    '''
    # --- 1. 验证输入并解包数据 ---
    required_keys = ['profile_data', 'y_coords', 'z_coords', 'lon_coords', 'lat_coords', 'projections', 'metadata']
    if not data_package or not all(k in data_package for k in required_keys):
        print(f"错误: 输入的 data_package 格式不完整。缺少键。")
        return
        
    profile_variable_2d = data_package['profile_data'].get(variable)
    if profile_variable_2d is None or np.all(profile_variable_2d.mask):
        print(f"警告: 变量 '{variable}' 的剖面数据无效，无法绘图。")
        return

    y_coords = data_package['y_coords']
    z_coords = data_package['z_coords']
    projections = data_package['projections']
    metadata = data_package['metadata']

    # --- 2. 在当前上下文中计算与DS相关的元数据 ---
    callers_local_vars = inspect.currentframe().f_back.f_locals
    ds_name_list = [var_name for var_name, var_val in callers_local_vars.items() if var_val is DS]
    if not ds_name_list:
        raise ValueError("无法在调用者环境中找到数据集变量名。请确保DS参数已正确传入。")
    ds_name = ds_name_list[0].upper()
    
    prop_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    eddy_color = prop_colors[1] if 'AC' in ds_name else prop_colors[0]

    # --- 3. 准备绘图元素 ---
    if profile_variable_2d.shape[1] != len(y_coords):
       y_coords = y_coords[:profile_variable_2d.shape[1]]

    Y_mesh, Z_mesh = np.meshgrid(y_coords, z_coords)
    
    if variable == 'vorticity': cbar_label, cmap, clim = r'$\zeta/f$', 'seismic', (-0.3, 0.3)
    elif variable in ['thetao']: cbar_label, cmap = 'Temperature (°C)', 'rainbow'
    elif variable in ['salinity', 'so']: cbar_label, cmap = 'Salinity (psu)', 'viridis'
    elif variable in ['u', 'v', 'uo', 'vo']: cbar_label, cmap = 'Velocity (m/s)', 'RdBu_r'
    else: cbar_label, cmap, clim = variable, 'viridis', None
    
    if 'clim' not in locals() or clim is None:
        valid_values = profile_variable_2d[~profile_variable_2d.mask]
        clim = (valid_values.min(), valid_values.max()) if valid_values.size > 0 else (0,1)
        if variable in ['u', 'v', 'uo', 'vo']:
            max_abs = np.max(np.abs(valid_values)) if valid_values.size > 0 else 1
            clim = (-max_abs, max_abs)

    # --- 4. 执行绘图 ---
    fig, ax = plt.subplots(figsize=(20, 15))
    
    k_val, b_val = metadata.get('k'), metadata.get('b')
    if k_val is not None and b_val is not None:
        title = (f"Vertical Profile of {cbar_label} for Track {ds_name}{metadata['eddy_no']} "
                 f"on {metadata['date_str']}, y={k_val:.2f}x+{b_val:.2f}")
    else:
         title = (f"Vertical Profile of {cbar_label} for Track {ds_name}{metadata['eddy_no']} "
                 f"on {metadata['date_str']}")

    ax.set_title(title, fontsize=20)
    ax.set_xlabel('Distance from Eddy Center Projection (km)', fontsize=18)
    ax.set_ylabel('Depth (m)', fontsize=18)
    ax.tick_params(labelsize=14)
    
    pc = ax.pcolormesh(Y_mesh, Z_mesh, profile_variable_2d, cmap=cmap, shading='auto', vmin=clim[0], vmax=clim[1])
    cbar = fig.colorbar(pc, ax=ax, orientation='vertical', fraction=0.046, pad=0.04)
    cbar.set_label(cbar_label, fontsize=18)
    cbar.ax.tick_params(labelsize=14)

    ax.axvline(0, color='black', linestyle='--', linewidth=2, label='Eddy Center Projection')
    for i, dist in enumerate(projections['radius']): 
        ax.axvline(dist, color='r', linestyle='--', linewidth=2, label='Effective Radius Projection' if i == 0 else "")
    for i, dist in enumerate(projections['contour']): 
        ax.axvline(dist, color=eddy_color, linestyle=':', linewidth=2, label='Effective Contour Projection' if i == 0 else "")

    ax.set_ylim(z_coords.max(), z_coords.min())
    if xmin is not None and xmax is not None: ax.set_xlim(xmin, xmax)
    if ymin is not None and ymax is not None: ax.set_ylim(ymax, ymin)
    
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), fontsize=18)

    # --- 5. 保存和显示 ---
    if save_fig:
        output_dir = "plot_vertical_glorys"
        os.makedirs(output_dir, exist_ok=True)
        date_fn = metadata['date_str'].replace('-', '')
        
        if k_val is not None and b_val is not None:
            k_str = f"k{k_val:.2f}"
            b_str = f"b{b_val:.2f}"
            base_filename = (f"{ds_name}{metadata['eddy_no']}_vertical_{variable}_{date_fn}_{k_str}{b_str}.png")
        else:
            base_filename = (f"{ds_name}{metadata['eddy_no']}_vertical_{variable}_{date_fn}.png")

        plt.savefig(os.path.join(output_dir, base_filename), dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {os.path.join(output_dir, base_filename)}")

    if show_fig:
        plt.show()

    plt.close(fig)

def get_raytraceR_inputs(data_package: dict,
                                f_coriolis: float, # 科里奥利频率 f (例如 1.454e-4)。**必须提供**。
                                rho0: float = 1025.0, # 参考密度 (kg/m^3)。
                                g: float = 9.81, # 重力加速度 (m/s^2)。
                                m0: float = 1/100.0, # 初始垂直波数。
                                omega_factor: float = 0.97, # 固有频率 omega = omega_factor * f。
                                v0_amplitude: float = 0.1, # 初始波速振幅。
                                thresh_val: float = 0.1e-2, # 内部反射阈值。
                                chstart_val: int = 1, # 初始特征线 (1或2)。
                                filter_sigma_z: float = 0.5, # 沿Z方向高斯平滑的标准差（用于导数计算，防止噪声放大）。
                                filter_sigma_y: float = 0.5 # 沿Y方向高斯平滑的标准差。
                               ) -> dict:
    """
    将来自 track.regrid_vertical_slice 的数据包转换为 MATLAB raytraceR 函数所需的输入格式。

    假设输入 data_package 的结构如下（来自 track.get_vertical_glorys 和 regrid_vertical_slice 处理后）：
    - z_coords (np.ndarray): 等间距的深度 (Z) 坐标，单位为米(m)。
    - y_coords (np.ndarray): 等间距的横流 (Y) 坐标，单位为公里(km)。
    - lon_coords (np.ndarray): 剖面线上每个点对应的经度。
    - lat_coords (np.ndarray): 剖面线上每个点对应的纬度。
    - profile_data (dict): 包含以下二维 np.ndarray (或 np.ma.MaskedArray) 数组的字典:
        - 'u': 纬向速度 (uo)
        - 'v': 经向速度 (vo)
        - 'salinity': 盐度 (Practical Salinity, SP)
        - 'thetao': 位势温度 (potential temperature, pt)
    
    参数:
        data_package (dict): 包含 GLORYS 垂直切片数据的字典，通常是
                             track.regrid_vertical_slice 的输出。
        f_coriolis (float): 科里奥利频率 f (例如 1.454e-4 rad/s)。**必须提供**。
        rho0 (float): 参考密度 (kg/m^3)。
        g (float): 重力加速度 (m/s^2)。
        m0 (float): 初始垂直波数。
        omega_factor (float): 固有频率的倍数，omega = omega_factor * f。
        v0_amplitude (float): 初始波速振幅。
        thresh_val (float): 内部反射阈值。
        chstart_val (int): 初始特征线 (1或2)。
        filter_sigma_z (float): 沿深度方向高斯平滑的标准差（用于导数计算）。
        filter_sigma_y (float): 沿横流方向高斯平滑的标准差。

    返回:
        dict: 包含 raytraceR 所需所有输入参数的字典。
    """

    # 1. 提取和映射输入数据
    # raytraceR 的 Z 是垂直/深度，Y 是水平横流。
    z_g = data_package['z_coords'] # 深度坐标，单位：米 (m)
    y_g_km = data_package['y_coords'] # 横流坐标，单位：公里 (km)
    lon_coords_1d = data_package['lon_coords'] # 剖面经度 (1D array)
    lat_coords_1d = data_package['lat_coords'] # 剖面纬度 (1D array)

    # 将横流坐标从公里转换为米，以保持 raytraceR 内部单位一致性
    y_g = y_g_km * 1000.0 # 横流坐标，单位：米 (m)

    # 从 profile_data 中获取核心物理量
    # 它们已经是 (深度点数, 横流点数) 的形状
    # .filled(np.nan) 用于处理 MaskedArray，将其掩码值填充为 NaN
    u_data = data_package['profile_data']['u'].filled(np.nan)
    v_data = data_package['profile_data']['v'].filled(np.nan)
    salinity_data = data_package['profile_data']['salinity'].filled(np.nan) # Practical Salinity (SP)
    thetao_data = data_package['profile_data']['thetao'].filled(np.nan) # Potential Temperature (pt)

    # 检查维度和数据有效性
    expected_shape_2d = (len(z_g), len(y_g_km))
    if not (u_data.shape == salinity_data.shape == thetao_data.shape == expected_shape_2d):
        raise ValueError(f"Profile data shapes are inconsistent or do not match expected ({len(z_g)}, {len(y_g_km)}).")
    if not np.isfinite(f_coriolis):
        raise ValueError("Coriolis frequency 'f_coriolis' must be a finite number.")

    # 2. 计算网格间距
    dz = np.mean(np.diff(z_g)) # regrid_vertical_slice 保证了等间距
    dy = np.mean(np.diff(y_g)) # 横流坐标已经转换为米，这里计算米为单位的间距

    # 3. 对物理量进行平滑处理 (减少噪声对导数计算的影响)
    # 确保 sigma 值与数据尺寸兼容
    u_smoothed = gaussian_filter(u_data, sigma=(filter_sigma_z, filter_sigma_y))
    v_smoothed = gaussian_filter(v_data, sigma=(filter_sigma_z, filter_sigma_y))
    salinity_smoothed = gaussian_filter(salinity_data, sigma=(filter_sigma_z, filter_sigma_y))
    thetao_smoothed = gaussian_filter(thetao_data, sigma=(filter_sigma_z, filter_sigma_y))
    
    # 处理平滑后可能出现的NaN（高斯滤波的边界效应），尤其是在mask区域外推的NaN
    u_smoothed[np.isnan(u_smoothed)] = 0.0 # 假设用0填充，或选择其他填充值
    v_smoothed[np.isnan(v_smoothed)] = 0.0
    salinity_smoothed[np.isnan(salinity_smoothed)] = np.nanmean(salinity_smoothed) if np.all(np.isnan(salinity_smoothed)) else salinity_smoothed[np.isnan(salinity_smoothed)].mean() # 用平均值填充，避免gsw出错
    thetao_smoothed[np.isnan(thetao_smoothed)] = np.nanmean(thetao_smoothed) if np.all(np.isnan(thetao_smoothed)) else thetao_smoothed[np.isnan(thetao_smoothed)].mean() # 同上

    # 4. 计算核心物理量: 密度, 浮力 (bg), N2, ug, F2, S2, s_M, omegamin
    
    # 4.1. 转换为 gsw 所需的 Absolute Salinity (SA) 和 Conservative Temperature (CT)
    # gsw 需要经度和纬度网格与数据形状匹配
    Z_mesh_m, Y_mesh_m = np.meshgrid(z_g, y_g, indexing='ij') # Z_mesh_m (深度, 米), Y_mesh_m (横流，米)
    
    # 将 1D lon/lat_coords 扩展为 2D 网格，匹配 (深度点数, 横流点数) 形状
    # 假设 lon_coords_1d 和 lat_coords_1d 是与 y_g_km 对应的，即沿横流方向的经纬度
    # 那么它们需要沿着深度方向广播
    lon_2d = np.tile(lon_coords_1d[np.newaxis, :], (len(z_g), 1))
    lat_2d = np.tile(lat_coords_1d[np.newaxis, :], (len(z_g), 1))

    # 计算压力（p_from_z 需要负深度）
    pressure_2d = gsw.p_from_z(-Z_mesh_m, lat_2d)

    # 实用盐度 -> 绝对盐度
    SA_processed = gsw.SA_from_SP(salinity_smoothed, pressure_2d, lon_2d, lat_2d)
    
    # 位温 -> 保守温度
    CT_processed = gsw.CT_from_pt(SA_processed, thetao_smoothed)

    # 计算位势密度 (referenced to 0 dbar)
    rho_potential_processed = gsw.rho(SA_processed, CT_processed, 0)
    
    # 计算浮力 bg = -g * (rho_potential - rho0) / rho0
    bg_processed = -g * (rho_potential_processed - rho0) / rho0

    # 4.2. 计算 N2 (浮力频率的平方)
    # gsw.Nsquared(SA, CT, p, lat=None, axis=0) 返回 (N2, p_mid)
    # axis=0 表示沿深度（第一个）维度计算梯度
    N2_output, p_mid_N2 = gsw.Nsquared(SA_processed, CT_processed, pressure_2d, lat=lat_2d, axis=0)
    
    # N2_output 的形状比原始数据在深度维度上少一个点 (因为是 mid-point 导数)
    # 需要将其插值或映射回原始深度网格的尺寸
    # 常用方法是插值回原始深度点或在第一个点进行外推/复制
    # 这里我们采用在第一个深度点复制 N2_output 的第一行，使其维度匹配
    N2_processed = np.vstack([N2_output[0:1, :], N2_output])
    # 再次裁剪以确保与原始深度维度完全一致（防止因浮点数问题多/少一行）
    N2_processed = N2_processed[:u_data.shape[0], :] 
    
    # 处理 N2 中可能出现的 NaN 或负值（物理上 N2 应为正或零）
    N2_processed[np.isnan(N2_processed)] = 0.0 # 用 0 填充 NaN
    N2_processed[N2_processed < 0] = 0.0 # 强制非负

    # 4.3. 计算 ug (背景地转流场)
    # 论文中 ug 是背景地转流的 x 分量。这里直接使用平滑后的 u 分量作为 ug。
    ug_processed = u_smoothed

    # 4.4. 计算 F2, S2 (根据 Whitt and Thomas 2013 论文定义)
    
    # zeta_g = -du_g/dy
    # np.gradient 自动处理 MaskedArray，但我们已填充为 NaN
    # 注意 dy 已经是米
    du_dy_smoothed = np.gradient(ug_processed, dy, axis=1) # du_g/dy
    zeta_g_processed = -du_dy_smoothed
    
    # F2 = f * (f + zeta_g)
    F2_processed = f_coriolis * (f_coriolis + zeta_g_processed)

    # S2 = f * du_g/dz
    # 注意 dz 已经是米
    du_dz_smoothed = np.gradient(ug_processed, dz, axis=0) # du_g/dz
    S2_processed = f_coriolis * du_dz_smoothed
    
    # 4.5. 计算 s_M = F2 / S2
    # 避免除以零或 NaN
    s_M_processed = np.divide(F2_processed, S2_processed, out=np.zeros_like(F2_processed), where=S2_processed!=0)
    s_M_processed[np.isnan(s_M_processed)] = 0.0 # 填充 NaN
    s_M_processed[np.isinf(s_M_processed)] = 0.0 # 填充 Inf

    # 4.6. 计算 omegamin
    # omegamin = sqrt(q/N2)
    # q = F2 * N2 - S2^4
    PV_q = F2_processed * N2_processed - S2_processed**4
    
    # 避免根号下出现负值或除以零
    omegamin_squared = np.divide(PV_q, N2_processed, out=np.zeros_like(PV_q), where=N2_processed!=0)
    omegamin_squared = np.maximum(0, omegamin_squared) # 确保非负

    omegamin_processed = np.sqrt(omegamin_squared)
    omegamin_processed = np.real(omegamin_processed) # 确保是实数
    omegamin_processed[np.isnan(omegamin_processed)] = 0.0 # 填充 NaN
    omegamin_processed[np.isinf(omegamin_processed)] = 0.0 # 填充 Inf

    # 5. 二阶导数 (根据论文示例，默认可以设为零，或者在这里计算)
    # 鉴于计算二阶导数对噪声敏感且需要额外代码，我们这里保持为零，
    # 除非用户明确要求并提供计算方法。
    d2udy2_processed = np.zeros_like(u_data)
    d2udz2_processed = np.zeros_like(u_data)
    d2bdy2_processed = np.zeros_like(u_data)
    d2bdz2_processed = np.zeros_like(u_data)
    d2bdzdy_processed = np.zeros_like(u_data)
    d2udzdy_processed = np.zeros_like(u_data)

    # 6. 计算领域长度和波参数
    omega = omega_factor * f_coriolis
    Lz = z_g.max() - z_g.min()
    Ly = y_g.max() - y_g.min() # y_g 已经是米

    # m0 = 1/minv，所以 minv = 1/m0
    minv_val = 1.0 / m0

    # 准备输出字典，包含 raytraceR 所需所有参数
    raytraceR_inputs = {
        'Z': -Z_mesh_m, # 2D 深度网格 (深度, 横流)，单位：米 (m)
        'dz': dz, # 深度间距 (米)
        'Y': Y_mesh_m, # 2D 横流网格 (深度, 横流)，单位：米 (m)
        'dy': dy, # 横流间距 (米)
        'F2': F2_processed, # 2D 矩阵 (深度, 横流)
        'S2': S2_processed, # 2D 矩阵 (深度, 横流)
        'N2': N2_processed, # 2D 矩阵 (深度, 横流)
        'ug': ug_processed, # 2D 矩阵 (深度, 横流)
        'bg': bg_processed, # 2D 矩阵 (深度, 横流)
        's_M': s_M_processed, # 2D 矩阵 (深度, 横流)
        'omegamin': omegamin_processed, # 2D 矩阵 (深度, 横流)
        'f': f_coriolis, # 标量
        'g': g, # 标量
        'rho0': rho0, # 标量
        'omega': omega, # 标量
        'minv': minv_val, # 标量
        'v0': v0_amplitude, # 标量
        
        # 二阶导数 (目前设为零)
        'd2udy2': d2udy2_processed,
        'd2udz2': d2udz2_processed,
        'd2bdy2': d2bdy2_processed,
        'd2bdz2': d2bdz2_processed,
        'd2bdzdy': d2bdzdy_processed,
        'd2udzdy': d2udzdy_processed,
        
        # 领域长度 (根据 Z, Y 范围计算)
        'Lz': Lz,
        'Ly': Ly,
        
        # 其他 raytraceR 需要但通常在调用时设定的参数
        'thresh': thresh_val,
        'm0': m0, # 与 minv 互补，此处也返回
        'chstart': chstart_val,
        # projections数据用于Matlab中绘图
        'eddy_radius': data_package['projections'].get('radius', []), # 半径投影
        'eddy_contour': data_package['projections'].get('contour', []), # 等值线投影
    }
    
    print("\n--- Data Preparation Summary for raytraceR ---")
    print(f"Z (Depth) meshgrid shape: {Z_mesh_m.shape}, min/max: {Z_mesh_m.min():.1f}/{Z_mesh_m.max():.1f} m")
    print(f"Y (Cross-stream) meshgrid shape: {Y_mesh_m.shape}, min/max: {Y_mesh_m.min()/1000:.1f}/{Y_mesh_m.max()/1000:.1f} km")
    print(f"dz (Depth spacing): {dz:.3f} m")
    print(f"dy (Cross-stream spacing): {dy:.3f} m")
    print(f"Calculated F2 shape: {F2_processed.shape}, min/max: {np.nanmin(F2_processed):.2e}/{np.nanmax(F2_processed):.2e}")
    print(f"Calculated S2 shape: {S2_processed.shape}, min/max: {np.nanmin(S2_processed):.2e}/{np.nanmax(S2_processed):.2e}")
    print(f"Calculated N2 shape: {N2_processed.shape}, min/max: {np.nanmin(N2_processed):.2e}/{np.nanmax(N2_processed):.2e}")
    print(f"Calculated ug shape: {ug_processed.shape}, min/max: {np.nanmin(ug_processed):.2e}/{np.nanmax(ug_processed):.2e}")
    print(f"Calculated bg shape: {bg_processed.shape}, min/max: {np.nanmin(bg_processed):.2e}/{np.nanmax(bg_processed):.2e}")
    print(f"Calculated omegamin shape: {omegamin_processed.shape}, min/max: {np.nanmin(omegamin_processed):.2e}/{np.nanmax(omegamin_processed):.2e}")
    print(f"Coriolis frequency (f): {f_coriolis:.4e} rad/s")
    print(f"Wave frequency (omega): {omega:.4e} rad/s")
    print("--------------------------------------------\n")

    return raytraceR_inputs

def init_worker(eddy_data_shared: dict):
    """
    (Initializer函数) 为multiprocessing的子进程初始化共享的、只读的数据。

    功能:
        这个函数在每个工作进程启动时仅被调用一次。它接收大的数据集
        (如涡旋数据字典)并将其设置为该进程的全局变量。这可以极大地避免
        在每个任务间重复传输大数据的开销，是性能优化的关键。
    """
    global worker_eddy_datasets
    worker_eddy_datasets = eddy_data_shared

def check_single_track(track_data, argo_points_by_date, start_date, end_date, ds_name):
    """
    (内部辅助函数) 检查单个涡旋轨迹是否与Argo数据有交集。

    功能:
        这是一个纯计算函数，会为所有在时间范围内的涡旋返回结果，
        并附带一个布尔标志来说明其是否与Argo浮标有交集。

    参数:
        track_data (list): 单条涡旋的轨迹数据 (list of lists)。
        argo_points_by_date (dict): 按日期组织的Argo位置字典。
        start_date (pd.Timestamp): 检查的开始日期。
        end_date (pd.Timestamp): 检查的结束日期。
        ds_name (str): 数据集名称 (如 'ACS')。

    返回:
        dict | None: 如果涡旋在时间范围内，返回包含绘图信息的字典，否则返回None。
    """
    num, time, center_lon, center_lat, _, _, contour_lon, contour_lat, radius, _, _ = zip(*track_data)
    dates = convert_date(time)
    
    indices_in_range = np.where((dates >= start_date) & (dates <= end_date))[0]
    if indices_in_range.size == 0:
        return None

    has_interaction = False
    contours_to_plot = []
    for i in indices_in_range:
        current_date = dates[i].normalize()
        if current_date in argo_points_by_date:
            contour_poly = Polygon(zip(contour_lon[i], contour_lat[i]))
            for argo_point in argo_points_by_date[current_date]:
                inside_poly = contour_poly.contains(argo_point)
                center = np.array([center_lon[i], center_lat[i]])
                point_coord = np.array([argo_point.x, argo_point.y])
                distance = np.linalg.norm(point_coord - center)
                inside_circle = distance <= (radius[i] / 111320) * circle_enlargement_factor
                if inside_poly or inside_circle:
                    contours_to_plot.append((contour_lon[i], contour_lat[i]))
                    has_interaction = True
                    break
    
    in_range_segments = []
    splits = np.where(np.diff(indices_in_range) != 1)[0] + 1
    contiguous_blocks = np.split(indices_in_range, splits)
    for block in contiguous_blocks:
        if block.size > 0:
            in_range_segments.append((np.array(center_lon)[block], np.array(center_lat)[block]))

    start_idx = indices_in_range[0]
    text_info = {"text": f"{ds_name}{num[0]}", "lon": center_lon[start_idx], "lat": center_lat[start_idx]}

    return {
        "ds_name": ds_name, "center_lon": center_lon, "center_lat": center_lat,
        "in_range_segments": in_range_segments, "contours_to_plot": contours_to_plot,
        "text_info": text_info, "is_ace": 'AC' in ds_name.upper(),
        "has_interaction": has_interaction
    }

def plot_all_tracks_in_range(
    start_date_str: str,
    end_date_str: str,
    eddy_datasets: dict | None = None,
    plot_unrelated_eddies: bool = False,
    plot_unrelated_argo: bool = True,
    save_fig: bool = False,
    show_fig: bool = True,
    do_threshold: float = 50.0,
    salinity_threshold: float = 0.0,
    temperature_threshold: float = 0.0,
    depth_interval: float = 100.0,
    depth_merge_tolerance: float = 10.0,
    duplicate_depth_strategy: str = 'best_qc',
    anomaly_min_depth: float | None = 300.0,
    anomaly_color_by: str = 'delta_do',
    fix_delta_do_colorbar: bool = True,
    delta_do_cbar_min: float = 50.0,
    delta_do_cbar_max: float = 100.0,
    delta_do_cbar_ticks: list | None = None
):
    """(核心绘图) 指定时间段内涡旋轨迹 + Argo ΔDO 异常代表点（仅采用 ΔDO 方法）。

    工作流程：
      1. 装载时间范围内 Argo 数据 → 过滤地理范围 → 计算 ΔDO 异常。
      2. 可选按 anomaly_min_depth 过滤异常深度。
      3. 每个剖面保留 delta_do（或 do_value）最大的一条。
      4. 按 anomaly_color_by 着色：'delta_do' (默认) 或 'do_value'。

    参数:
        start_date_str, end_date_str: 日期范围。
        eddy_datasets: 预加载涡旋数据；并行 worker 模式下从全局读取。
        plot_unrelated_eddies: 是否绘制未与 Argo 交互的涡旋。
        plot_unrelated_argo: 是否额外绘制所有 Argo 剖面位置（空心圆），用于提供基准分布背景。
        save_fig, show_fig: 输出控制。
        do_threshold / salinity_threshold / temperature_threshold / depth_interval / depth_merge_tolerance / duplicate_depth_strategy: 传给 calculate_delta_do。
        anomaly_min_depth: (可选) 仅保留异常深度 >= 此值；None 不限制。
        anomaly_color_by: 'delta_do' 或 'do_value'。
        fix_delta_do_colorbar: 若为 True 且按 delta_do 着色，则强制使用 [delta_do_cbar_min, delta_do_cbar_max] 作为色标范围。
        delta_do_cbar_min / delta_do_cbar_max: ΔDO 色标固定范围上下限（仅在 fix_delta_do_colorbar=True 且 anomaly_color_by='delta_do' 时生效）。
        delta_do_cbar_ticks: 自定义 ΔDO 色标刻度列表（None 自动：若只提供上下限则显示两端；若范围>30 添加中点）。
    """
    # --- 0. 确定数据源 ---
    local_eddy_datasets = eddy_datasets
    # 如果作为并行worker运行，eddy_datasets会是None，此时从全局变量获取
    if local_eddy_datasets is None:
        global worker_eddy_datasets
        local_eddy_datasets = worker_eddy_datasets

    start_date = pd.to_datetime(start_date_str)
    end_date = pd.to_datetime(end_date_str)
    
    # --- 1. 动态加载Argo数据 ---
    min_year, max_year = start_date.year, end_date.year
    argo_data_list = []
    for year in range(min_year, max_year + 1):
        try:
            argo_data_list.append(load_argo_data(year, data_dir=argo_path))
        except FileNotFoundError:
            print(f"Warning: Argo data for year {year} not found, skipping.")
    if not argo_data_list:
        print(f"No Argo data found for the period {start_date_str} to {end_date_str}.")
        return

    argo_data = pd.concat(argo_data_list, ignore_index=True)
    
    # --- 2. 预处理Argo数据 ---
    argo_dates = pd.to_datetime(argo_data[['Year', 'Month', 'Day']])
    argo_in_range = argo_data[(argo_dates >= start_date) & (argo_dates <= end_date)]
    needed_argo_data = pd.DataFrame()
    base_argo_positions = pd.DataFrame()
    if not argo_in_range.empty:
        geo_mask = (
            (argo_in_range['Longitude'] >= lonmin) & (argo_in_range['Longitude'] <= lonmax) &
            (argo_in_range['Latitude'] >= latmin) & (argo_in_range['Latitude'] <= latmax)
        )
        argo_in_geo_range = argo_in_range[geo_mask].copy()
        if not argo_in_geo_range.empty:
            # 所有剖面基础位置（避免按深度重复）
            base_argo_positions = (
                argo_in_geo_range.sort_values(['Profile_number','Depth'])
                .groupby('Profile_number', as_index=False)
                .first()[['Profile_number','Longitude','Latitude']]
            )
            anomalies = calculate_delta_do(
                argo_in_geo_range,
                depth_interval=depth_interval,
                do_threshold=do_threshold,
                salinity_threshold=salinity_threshold,
                temperature_threshold=temperature_threshold,
                depth_merge_tolerance=depth_merge_tolerance,
                duplicate_depth_strategy=duplicate_depth_strategy,
                remove_outliers=True,
                verbose=False
            )
            if not anomalies.empty:
                if anomaly_min_depth is not None:
                    anomalies = anomalies[anomalies['depth'] >= anomaly_min_depth]
                if not anomalies.empty:
                    sort_field = 'delta_do' if (
                        anomaly_color_by == 'delta_do' and 'delta_do' in anomalies.columns
                    ) else 'do_value'
                    anomalies_sorted = anomalies.sort_values(by=[sort_field], ascending=False)
                    anomalies_unique = anomalies_sorted.drop_duplicates(subset='Profile_number', keep='first')
                    needed_argo_data = anomalies_unique.rename(columns={'depth': 'Anomaly_depth'})

    argo_points_by_date = defaultdict(list)
    if not needed_argo_data.empty:
        for _, row in needed_argo_data.iterrows():
            date_key = pd.Timestamp(year=int(row['Year']), month=int(row['Month']), day=int(row['Day']))
            argo_points_by_date[date_key].append(Point(row['Longitude'], row['Latitude']))

    # --- 3. 检查所有涡旋轨迹 ---
    all_tracks_with_names = [(track, ds_name) for ds_name, ds_data in local_eddy_datasets.items() for track in ds_data]
    results = [check_single_track(track, argo_points_by_date, start_date, end_date, ds_name) for track, ds_name in all_tracks_with_names]

    # --- 4. 绘图 ---
    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    fig, ax = plt.subplots(figsize=(40, 30))
    ax.set_title(
        f"Eddy Tracks and Argo ΔDO Anomalies ({start_date.strftime('%Y-%m-%d')} to {end_date.strftime('%Y-%m-%d')})",
        fontsize=20
    )
    ax.set_xlabel('Longitude', fontsize=20); ax.set_ylabel('Latitude', fontsize=20)
    world.plot(color='lightgrey', edgecolor='white', ax=ax)
    
    prop_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    eddy_colors = {'ACE': prop_colors[1], 'CE': prop_colors[0]}

    for result in filter(None, results):
        has_interaction = result['has_interaction']
        if not has_interaction and not plot_unrelated_eddies:
            continue
        
        is_ace = result['is_ace']
        color = eddy_colors['ACE'] if is_ace else eddy_colors['CE']
        text_color = 'red' if has_interaction else 'black'
        
        ax.plot(result['center_lon'], result['center_lat'], color=color, alpha=0.4, linestyle='--', zorder=4)
        for lon_seg, lat_seg in result['in_range_segments']:
            ax.plot(lon_seg, lat_seg, color=color, alpha=0.8, linestyle='-', zorder=4)
        
        info = result['text_info']
        ax.text(info['lon'], info['lat'], info['text'], fontsize=12, color=text_color, weight='bold', zorder=5)

        if has_interaction:
            for contour_lon, contour_lat in result['contours_to_plot']:
                ax.plot(contour_lon, contour_lat, color=color, linewidth=1, alpha=0.5, zorder=4, linestyle=':')

    if plot_unrelated_argo and not base_argo_positions.empty:
        ax.scatter(
            base_argo_positions['Longitude'], base_argo_positions['Latitude'],
            facecolors='none', edgecolors='gray', linewidths=0.8, s=36,
            label='All Argo Profiles (baseline)', zorder=2
        )

    if not needed_argo_data.empty:
        if anomaly_color_by == 'delta_do' and 'delta_do' in needed_argo_data.columns:
            # 应用固定色标（可选）
            scatter_kwargs = {}
            if fix_delta_do_colorbar:
                scatter_kwargs.update(dict(vmin=delta_do_cbar_min, vmax=delta_do_cbar_max))
            sc = ax.scatter(
                needed_argo_data['Longitude'], needed_argo_data['Latitude'],
                c=needed_argo_data['delta_do'], cmap='Reds', s=70,
                edgecolors='black', linewidths=0.5,
                label=f'ΔDO ≥ {do_threshold} μmol kg⁻¹ @ depth ≥ {anomaly_min_depth} m', zorder=3,
                **scatter_kwargs
            )
            cbar = plt.colorbar(sc, ax=ax, orientation='horizontal', fraction=0.046, pad=0.04)
            cbar.set_label('ΔDO / μmol·kg⁻¹', fontsize=20); cbar.ax.tick_params(labelsize=14)
            # 自动或用户自定义刻度
            if fix_delta_do_colorbar:
                if delta_do_cbar_ticks is not None:
                    cbar.set_ticks(delta_do_cbar_ticks)
                else:
                    # 默认：只显示上下限；若范围较大则加中点
                    rng = delta_do_cbar_max - delta_do_cbar_min
                    if rng > 30:
                        mid = (delta_do_cbar_max + delta_do_cbar_min) / 2
                        cbar.set_ticks([delta_do_cbar_min, mid, delta_do_cbar_max])
                    else:
                        cbar.set_ticks([delta_do_cbar_min, delta_do_cbar_max])
        else:
            color_values = needed_argo_data.get('do_value') if 'do_value' in needed_argo_data.columns else needed_argo_data.get('DO')
            sc = ax.scatter(
                needed_argo_data['Longitude'], needed_argo_data['Latitude'],
                c=color_values, cmap='bwr', s=60, vmin=150, vmax=240,
                edgecolors='black', linewidths=0.5,
                label='Argo DO Anomaly Profiles',
                zorder=3
            )
            cbar = plt.colorbar(sc, ax=ax, orientation='horizontal', fraction=0.046, pad=0.04)
            cbar.set_label('DO / μmol·kg⁻¹', fontsize=20); cbar.ax.tick_params(labelsize=14)

    ax.set_xlim(lonmin, lonmax); ax.set_ylim(latmin, latmax)
    ax.tick_params(axis='both', which='major', labelsize=16); ax.set_aspect('equal')
    
    legend_elements = [
        Line2D([0], [0], color=eddy_colors['ACE'], lw=2, label='ACE Track'),
        Line2D([0], [0], color=eddy_colors['CE'], lw=2, label='CE Track')
    ]
    handles, labels = ax.get_legend_handles_labels()
    extra_labels = [f"Argo ΔDO ≥ {do_threshold}", 'Argo DO Anomaly Profiles', 'All Argo Profiles (baseline)']
    added = None
    for lab in extra_labels:
        if lab in labels:
            added = handles[labels.index(lab)]
            break
    if added is not None:
        ax.legend(handles=legend_elements + [added], fontsize=18, loc='upper left')
    else:
        ax.legend(handles=legend_elements, fontsize=18, loc='upper left')

    # --- 5. 输出控制 ---
    if save_fig:
        output_dir = Path("plot_all_tracks_in_range")
        output_dir.mkdir(exist_ok=True, parents=True)
        base_filename = f"All_Tracks_{start_date_str}_to_{end_date_str}.png"
        save_path = output_dir / base_filename
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {save_path}")
    
    if show_fig:
        plt.show()
    
    plt.close(fig)

    # 返回本时间段内发生交互(标红)的涡旋标签列表（如 ACxxx / CEyyy）
    interacted_eddies = []
    for r in filter(None, results):
        if r['has_interaction']:
            interacted_eddies.append(r['text_info']['text'])
    return interacted_eddies

def worker_wrapper(args: tuple):
    """multiprocessing worker 包装函数。args: (start_date_str, end_date_str, plot_unrelated_eddies). 返回本月交互涡旋标签列表。"""
    start_d, end_d, unrelated_flag = args
    try:
        return plot_all_tracks_in_range(
            start_date_str=start_d,
            end_date_str=end_d,
            plot_unrelated_eddies=unrelated_flag,
            save_fig=True,
            show_fig=False,
        )
    except Exception as e:
        print(f"!!! ERROR processing period {start_d}: {e}")
        return []

def run_batch_plotting_multiprocessing(
    start_date_str: str,
    end_date_str: str,
    eddy_datasets: dict,
    num_workers: int,
    plot_unrelated_eddies: bool = False
):
    """
    (批处理控制器) 使用multiprocessing启动一个扁平化的并行绘图作业，并带有进度条。
    
    功能:
        1. 创建一个横跨指定数量核心的进程池。
        2. 使用initializer高效地将涡旋数据共享给所有工作进程。
        3. 按月份切分任务，并使用tqdm实时显示处理进度。

    参数:
        start_date_str (str): 批处理的开始日期 'YYYY-MM-DD'。
        end_date_str (str): 批处理的结束日期 'YYYY-MM-DD'。
        eddy_datasets (dict): 【已加载】的、将被共享给所有进程的涡旋数据集。
        num_workers (int): 需要启动的并行工作进程数（核心数）。
        plot_unrelated_eddies (bool): 是否在批处理中绘制无关涡旋。
    """
    print("="*60)
    print("      Multiprocessing Batch Plotting with Progress Bar      ")
    print("="*60)
    print(f"Strategy: Creating a single pool of {num_workers} workers.")
    
    # --- 1. 创建按月切分的任务列表 ---
    month_starts = pd.date_range(start=start_date_str, end=end_date_str, freq='MS')
    tasks = [
        (
            start_date.strftime('%Y-%m-%d'), 
            (start_date + pd.tseries.offsets.MonthEnd(1)).strftime('%Y-%m-%d'),
            plot_unrelated_eddies
        ) 
        for start_date in month_starts
    ]
    
    print(f"[*] Created {len(tasks)} monthly plotting tasks to be processed by {num_workers} cores.")
    
    # --- 2. 启动进程池并执行任务 ---
    start_time_total = tm.time()
    
    # 使用 initializer 来高效地传递一次大的涡旋数据
    collected = []
    with multiprocessing.Pool(processes=num_workers, initializer=init_worker, initargs=(eddy_datasets,)) as pool:
        print("[*] Processing tasks...")
        for month_result in tqdm(pool.imap_unordered(worker_wrapper, tasks), total=len(tasks)):
            if month_result:
                collected.extend(month_result)
        
    end_time_total = tm.time()
    total_duration_minutes = (end_time_total - start_time_total) / 60
    
    unique_interacted = sorted(set(collected))
    print("\n" + "="*60)
    print("--- All Plotting Tasks Have Finished ---")
    print(f"Total execution time: {total_duration_minutes:.2f} minutes.")
    print(f"Total interacted eddies: {len(unique_interacted)}")
    if unique_interacted:
        preview = ", ".join(unique_interacted[:20])
        print(f"Sample (first 20): {preview}{' ...' if len(unique_interacted)>20 else ''}")
    print("="*60)
    return unique_interacted

def plot_argo_hotspots(
    start_year: int,
    end_year: int,
    depth_threshold: float = 500.0,
    save_fig: bool = False,
    show_fig: bool = True
):
    """
    绘制指定时间范围内所有Argo数据中，满足特定条件的“热点”分布图。

    此函数旨在识别并可视化高溶解氧（DO）的区域。它会找出每个Argo剖面在指定深度
    之下的最大溶解氧点，并将其绘制在地图上。

    参数:
        start_year (int): 开始的年份。
        end_year (int): 结束的年份。
        depth_threshold (float, optional): 筛选的最小深度 (Depth >= threshold)。默认为 500.0。
        save_fig (bool, optional): 是否将生成的图像保存到文件。默认为 False。
        show_fig (bool, optional): 是否在屏幕上显示生成的图像。默认为 True。
    """
    print("--- Starting Deep Argo DO Hotspot Analysis ---")
    
    # --- 1. 循环加载并合并所有年份的数据 ---
    all_yearly_data = []
    print(f"Loading all Argo data from {start_year} to {end_year}...")
    for year in range(start_year, end_year + 1):
        print(f"Processing year: {year}...")
        try:
            yearly_argo_data = load_argo_data(year=year)
            all_yearly_data.append(yearly_argo_data)
        except FileNotFoundError as e:
            print(f"Warning: {e}. Skipping year {year}.")
        except Exception as e:
            print(f"An error occurred processing year {year}: {e}")

    if not all_yearly_data:
        print("No Argo data was loaded. Aborting plot.")
        return

    # --- 2. 核心筛选逻辑：先筛选深度，再找DO最高点 ---
    print("Combining all data and applying filters...")
    combined_df = pd.concat(all_yearly_data, ignore_index=True)
    
    # 步骤 A: 首先，只保留深度大于等于阈值的数据
    deep_argo_df = combined_df[combined_df['Depth'] >= depth_threshold].copy()
    
    if deep_argo_df.empty:
        print(f"No Argo data found below {depth_threshold}m. Aborting plot.")
        return

    # 步骤 B: 然后，在这个深层数据子集上，找出每个浮标DO最高的点
    print("Finding the highest DO record below threshold for each profile...")
    
    # --- 核心修正点 ---
    # 1. 先计算出每组最大DO值的行索引
    max_do_indices = deep_argo_df.groupby('Profile_number')['DO'].idxmax()
    # 2. 清理索引：使用 .dropna() 去掉因整组都是NaN而产生的缺失值
    valid_indices = max_do_indices.dropna()

    if valid_indices.empty:
        print("No valid max DO points found after filtering. Aborting plot.")
        return
        
    # 3. 使用清理后的有效索引来选取行
    hotspot_points_df = deep_argo_df.loc[valid_indices].copy()
    # --- 修正结束 ---
        
    # --- 3. 准备绘图数据 ---
    print(f"Sorting {len(hotspot_points_df)} unique hotspot points for plotting...")
    hotspot_points_df.sort_values(by='DO', ascending=True, inplace=True)
    
    # --- 4. 开始绘图 ---
    print("Generating plot...")
    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    fig, ax = plt.subplots(figsize=(40, 30))
    ax.set_title(f'Highest DO Record of Each Argo Float (below {depth_threshold}m) from {start_year}-{end_year}', fontsize=20)
    ax.set_xlabel('Longitude', fontsize=20)
    ax.set_ylabel('Latitude', fontsize=20)
    world.plot(color='lightgrey', edgecolor='white', ax=ax)
    
    sc = ax.scatter(
        hotspot_points_df['Longitude'], 
        hotspot_points_df['Latitude'], 
        c=hotspot_points_df['DO'], 
        cmap='bwr',
        vmin=150,
        vmax=240,
        s=25,
        alpha=0.8,
        label=f'Highest DO point of each float (below {depth_threshold}m)'
    )
    
    cbar = plt.colorbar(sc, ax=ax, orientation='horizontal', fraction=0.046, pad=0.04)
    cbar.set_label('DO / μmol·kg⁻¹ at hotspot depth', fontsize=20)
    cbar.ax.tick_params(labelsize=14)

    # --- 5. 最终化绘图 ---
    ax.set_xlim(lonmin, lonmax)
    ax.set_ylim(latmin, latmax)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_aspect('equal')
    ax.legend(fontsize=18, loc='upper left')

    # --- 6. 保存和显示图片 ---
    if save_fig:
        output_dir = "argo_hotspot_plots"
        os.makedirs(output_dir, exist_ok=True)
        base_filename = f"Argo_DO_Hotspots_{start_year}_to_{end_year}.png"
        save_path = os.path.join(output_dir, base_filename)
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {save_path}")

    if show_fig:
        plt.show()

    plt.close(fig)
    print("--- Hotspot analysis finished. ---")

def calculate_delta_do(
    data: pd.DataFrame,
    depth_col: str = 'Depth',
    do_col: str = 'DO',
    salinity_col: str = 'Salinity',
    temperature_col: str = 'Temperature',
    depth_interval: float = 100.0,
    do_threshold: float = 50.0,
    salinity_threshold: float = 0.0,
    temperature_threshold: float = 0.0,
    depth_merge_tolerance: float = 10.0,
    duplicate_depth_strategy: str = 'best_qc',
    remove_outliers: bool = True,
    verbose: bool = False
) -> pd.DataFrame:
    """
    参考论文方法计算亚表层异常信号（以 DO 为主）。

    步骤概述：
    1. 按 Profile_number 分组，分别处理每个剖面；
    2. 计算 DO 随深度的一阶导数，定位 DO 的正峰值（代表可能的表层富氧水体俯冲信号）；
    3. 以 DO 峰深度为中心，取窗口 [p-Δp, p+Δp]，用两端点连线构造参考剖面；
    4. 在峰值同一深度计算 ΔDO、ΔSalinity 与 ΔTemperature（原始值减参考线值）；
    5. 以 ΔDO ≥ do_threshold 作为必要条件；如设置了 salinity_threshold 或 temperature_threshold > 0，可附加 |ΔSalinity/ΔTemperature| 过滤；
    6. 同一剖面内若有相距很近的多个候选深度（常见于峰值上下各一点），按 depth_merge_tolerance（dbar）合并，仅保留 delta_do 较大的记录。

    参数：
        data (pd.DataFrame): 包含多个剖面数据的表；需包含 Profile_number、深度、DO、盐度等列。
        depth_col (str): 深度列名，默认 'Depth'。
        do_col (str): 溶解氧列名，默认 'DO'。
        salinity_col (str): 盐度列名，默认 'Salinity'。
        temperature_col (str): 温度列名，默认 'Temperature'。
        depth_interval (float): 深度窗口半宽（dbar），默认 100.0（总宽度 2*Δp）。
        do_threshold (float): ΔDO 判定阈值（μmol/kg），默认 50.0。
        salinity_threshold (float): ΔSalinity 可选阈值（psu）；≤0 表示不启用盐度过滤，默认 0.0。
        temperature_threshold (float): ΔTemperature 可选阈值（°C）；≤0 表示不启用温度过滤，默认 0.0。
        depth_merge_tolerance (float): 同一 Profile 内“深度近邻合并”阈值（dbar）。若两个候选点深度差小于该值，仅保留 delta_do 较大的记录；设为 ≤0 表示不合并，默认 10.0。
        duplicate_depth_strategy (str): 处理同一剖面内“同深度多条记录”的策略。
            可选 'best_qc'|'first'|'mean'|'max'|'min'（默认 'best_qc'：按 DO 的 QC 优先级 1>2>5>8 选最佳；并列取首个）。
        remove_outliers (bool): 基础 QC 与规则过滤，默认 True。
        verbose (bool): 是否打印进度信息（开始处理/已处理/未检测到/总共检测到），默认 False。

    返回：
        pd.DataFrame: 每个满足条件的峰值一行，含
        Profile_number, depth, delta_do, delta_salinity, delta_temperature,
        do_value, salinity_value, temperature_value，
        以及 Year/Month/Day/Longitude/Latitude/Platform_number（若存在）。
        若无满足条件记录，返回空表。
    """
    
    # 检查必要的列是否存在
    required_cols = [depth_col, do_col, salinity_col, 'Profile_number']
    missing_cols = [col for col in required_cols if col not in data.columns]
    if missing_cols:
        print(f"警告：缺少必要的列: {missing_cols}")
        return pd.DataFrame()
    
    # 复制数据以避免修改原始数据
    input_data = data.copy()
    
    # 存储所有剖面的结果
    all_results = []
    
    # 按Profile_number分组处理
    profile_groups = input_data.groupby('Profile_number')
    total_profiles = len(profile_groups)
    if verbose:
        print(f"开始处理 {total_profiles} 个剖面...")
    
    processed_profiles = 0
    
    for profile_num, profile_data in profile_groups:
        # 质量控制：移除异常值和质量标记不良的数据
        if remove_outliers:
            # 应用Argo QC标准：仅保留等级为{1,2,5,8}的观测
            for var in [do_col, salinity_col, temperature_col]:
                qc_column_name = f"{var}_Flag"
                if qc_column_name in profile_data.columns:
                    good_qc_flags = ['1', '2', '5', '8', 1, 2, 5, 8]
                    bad_qc_mask = ~profile_data[qc_column_name].isin(good_qc_flags)
                    profile_data.loc[bad_qc_mask, var] = np.nan
            
            # 规则法：移除已知的错误值
            if do_col in profile_data.columns:
                bad_do_mask = profile_data[do_col] <= 1.0
                profile_data.loc[bad_do_mask, do_col] = np.nan
        
        # 移除包含NaN值的行
        drop_subset = [depth_col, do_col, salinity_col, temperature_col]
        profile_data_clean = profile_data.dropna(subset=drop_subset)
        
        if len(profile_data_clean) < 5:  # 需要至少5个数据点来计算导数和峰值
            continue  # 跳过数据点太少的剖面

        # 按深度排序
        profile_data_clean = profile_data_clean.sort_values(by=depth_col).reset_index(drop=True)

        # 处理重复或非严格递增的深度值以避免 np.gradient 内部出现除零 (dx1 或 dx2 = 0)
        # 策略：对同一深度的多条记录按 duplicate_depth_strategy 聚合为单条
        def _best_qc_pick(group: pd.DataFrame) -> pd.Series:
            qccol = f"{do_col}_Flag"
            priority = {1: 0, 2: 1, 5: 2, 8: 3}
            def rank(v):
                try:
                    iv = int(v)
                except Exception:
                    return 999
                return priority.get(iv, 999)
            if qccol in group.columns:
                ranks = group[qccol].apply(rank)
                min_rank = ranks.min()
                picked = group.loc[ranks[ranks == min_rank].index]
                # 并列时优先 DO 非空，随后保留第一条
                picked = picked[pd.notna(picked[do_col])] if do_col in picked.columns else picked
                return picked.iloc[0]
            # 无 QC 列则退化为 first
            return group.iloc[0]

        def _mean_pick(group: pd.DataFrame) -> pd.Series:
            # 对核心变量取均值，其他元数据取首个
            first = group.iloc[0].copy()
            for c in [do_col, salinity_col, temperature_col]:
                if c in group.columns:
                    first[c] = group[c].astype(float).mean()
            return first

        if profile_data_clean[depth_col].duplicated().any():
            strategy = (duplicate_depth_strategy or 'best_qc').lower()
            if strategy not in {'best_qc','first','mean','max','min'}:
                strategy = 'best_qc'
            grouped = list(profile_data_clean.groupby(depth_col, sort=False))
            picked_rows = []
            for depth_val, grp in grouped:
                if len(grp) == 1:
                    picked_rows.append(grp.iloc[0])
                    continue
                if strategy == 'best_qc':
                    picked_rows.append(_best_qc_pick(grp))
                elif strategy == 'first':
                    picked_rows.append(grp.iloc[0])
                elif strategy == 'mean':
                    picked_rows.append(_mean_pick(grp))
                elif strategy == 'max':
                    picked_rows.append(grp.loc[grp[do_col].idxmax()])
                elif strategy == 'min':
                    picked_rows.append(grp.loc[grp[do_col].idxmin()])
            profile_data_clean = pd.DataFrame(picked_rows)
            # 可能破坏原索引，重排并按深度排序
            profile_data_clean = profile_data_clean.sort_values(by=depth_col).reset_index(drop=True)
        depth_series = profile_data_clean[depth_col].values
        # 若仍非严格递增（可能存在逆序或噪声导致深度减小），尝试再次排序并去除 <= 前一值 的点
        # （二次排序保证顺序，随后用 np.diff > 0 过滤）
        if np.any(np.diff(depth_series) <= 0):
            # 先强制排序（已经排序过, 这里是保险）
            profile_data_clean = profile_data_clean.sort_values(by=depth_col).reset_index(drop=True)
            # 过滤掉与前一个深度差 <= 0 的观测
            cleaned_rows = [0]
            last_depth = profile_data_clean.loc[0, depth_col]
            for ridx in range(1, len(profile_data_clean)):
                dval = profile_data_clean.loc[ridx, depth_col]
                if dval > last_depth:  # 严格递增才保留
                    cleaned_rows.append(ridx)
                    last_depth = dval
            profile_data_clean = profile_data_clean.iloc[cleaned_rows].reset_index(drop=True)
            depth_series = profile_data_clean[depth_col].values
        # 若清理后点数不足或仍不严格递增则跳过
        if len(profile_data_clean) < 5 or (len(depth_series) > 1 and np.any(np.diff(depth_series) <= 0)):
            continue

        # 计算深度对DO和Salinity的斜率（一阶导数）
        depth_values = profile_data_clean[depth_col].values  # 已保证严格递增
        do_values = profile_data_clean[do_col].values
        salinity_values = profile_data_clean[salinity_col].values
        temperature_values = profile_data_clean[temperature_col].values

        # 使用中心差分法计算导数
        do_slopes = np.gradient(do_values, depth_values)

        # 重要说明：坐标系差异
        # 图像坐标系（海洋学习惯）：深度向下为正，图像上正斜率表示随深度增加变量增大
        # 数据计算坐标系：np.gradient计算dVar/dDepth，图像正斜率对应数值负斜率(<0)
        # 峰值识别：图像上斜率从正变负为负峰值，从负变正为正峰值
        # 因此代码中的判断条件与图像描述看似相反，但逻辑正确

        # 仅定位 DO 的峰值（保留正峰值作为俯冲信号候选）
        do_peaks = []
        for i in range(1, len(do_slopes) - 1):
            # DO 正峰值：图像上斜率从负变正（数值上从正变负）
            if do_slopes[i-1] > 0 and do_slopes[i+1] < 0:
                do_peaks.append((i, 'positive', depth_values[i]))
            # DO 负峰值：如需扩展也可纳入；当前仅保留正峰值

        if not do_peaks:
            continue

        # 对每个 DO 正峰，按同一深度计算 ΔDO、ΔSalinity、ΔTemperature
        use_salinity_filter = (salinity_threshold is not None and salinity_threshold > 0)
        use_temperature_filter = (temperature_threshold is not None and temperature_threshold > 0)
        profile_results = []

        for do_idx, do_type, target_depth in do_peaks:
            # 定义窗口 [p-Δp, p+Δp]
            depth_lower = max(0, target_depth - depth_interval)
            depth_upper = target_depth + depth_interval

            depth_mask = (depth_values >= depth_lower) & (depth_values <= depth_upper)
            if not np.any(depth_mask):
                continue

            depth_range = depth_values[depth_mask]
            do_range = do_values[depth_mask]
            salinity_range = salinity_values[depth_mask]
            temperature_range = temperature_values[depth_mask]

            if len(depth_range) < 2:
                continue

            # 线性参考剖面（两端点连线）
            do_ref_values = np.interp(
                depth_range, [depth_range[0], depth_range[-1]], [do_range[0], do_range[-1]]
            )
            salinity_ref_values = np.interp(
                depth_range, [depth_range[0], depth_range[-1]], [salinity_range[0], salinity_range[-1]]
            )
            temperature_ref_values = np.interp(
                depth_range, [depth_range[0], depth_range[-1]], [temperature_range[0], temperature_range[-1]]
            )

            # 取最接近目标深度的观测点
            target_idx = np.argmin(np.abs(depth_range - target_depth))

            delta_do = do_range[target_idx] - do_ref_values[target_idx]
            delta_salinity = salinity_range[target_idx] - salinity_ref_values[target_idx]
            delta_temperature = temperature_range[target_idx] - temperature_ref_values[target_idx]

            # 判定：ΔDO 为必需条件；ΔSalinity/ΔTemperature 在各自阈值>0时作为附加过滤（与条件）
            cond = (delta_do >= do_threshold)
            if use_salinity_filter:
                cond = cond and (not np.isnan(delta_salinity)) and (abs(delta_salinity) >= salinity_threshold)
            if use_temperature_filter:
                cond = cond and (not np.isnan(delta_temperature)) and (abs(delta_temperature) >= temperature_threshold)
            if cond:
                result = {
                    'Profile_number': profile_num,
                    'depth': target_depth,
                    'delta_do': delta_do,
                    'delta_salinity': delta_salinity,
                    'delta_temperature': delta_temperature,
                    'do_value': do_values[do_idx],
                    'salinity_value': salinity_range[target_idx],
                    'temperature_value': temperature_range[target_idx]
                }

                # 添加额外的剖面信息（如果存在）
                if 'Year' in profile_data_clean.columns:
                    result['Year'] = profile_data_clean['Year'].iloc[0]
                if 'Month' in profile_data_clean.columns:
                    result['Month'] = profile_data_clean['Month'].iloc[0]
                if 'Day' in profile_data_clean.columns:
                    result['Day'] = profile_data_clean['Day'].iloc[0]
                if 'Longitude' in profile_data_clean.columns:
                    result['Longitude'] = profile_data_clean['Longitude'].iloc[0]
                if 'Latitude' in profile_data_clean.columns:
                    result['Latitude'] = profile_data_clean['Latitude'].iloc[0]
                if 'Platform_number' in profile_data_clean.columns:
                    result['Platform_number'] = profile_data_clean['Platform_number'].iloc[0]

                profile_results.append(result)

        # 剖面内“深度近邻合并”：按 delta_do 降序贪心选取，避免在峰上下方重复取点
        if profile_results:
            if depth_merge_tolerance is not None and depth_merge_tolerance > 0:
                profile_results.sort(key=lambda r: (np.nan_to_num(r['delta_do'], nan=-np.inf)), reverse=True)
                kept = []
                kept_depths = []
                for rec in profile_results:
                    d = rec['depth']
                    if all(abs(d - kd) >= depth_merge_tolerance for kd in kept_depths):
                        kept.append(rec)
                        kept_depths.append(d)
                # 输出前按深度升序，便于查看
                kept.sort(key=lambda r: r['depth'])
                all_results.extend(kept)
            else:
                all_results.extend(profile_results)
        
        processed_profiles += 1
        if processed_profiles % 100 == 0 and verbose:
            print(f"已处理 {processed_profiles}/{total_profiles} 个剖面...")
    
    if not all_results:
        if verbose:
            print("未检测到满足阈值条件的DO异常信号。")
        return pd.DataFrame()
    
    results_df = pd.DataFrame(all_results)
    if verbose:
        print(f"总共检测到 {len(results_df)} 个潜在的DO异常信号，来自 {len(results_df['Profile_number'].unique())} 个剖面")

    return results_df