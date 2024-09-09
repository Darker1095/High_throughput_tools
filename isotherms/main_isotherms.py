import configparser
import math
import os
import re
import shutil
import threading
import time
from queue import Queue
from threading import Lock


class RASPA_Output_Data():
    '''
        RASPA输出文件对象
    '''
    '''
    示例：
        with open('./output.data','r') as f:
            str = f.read()
        output = RASPA_Output_Data(str)
        print(output.is_finished())
        print(output.get_absolute_adsorption())

    '''

    def __init__(self, output_string):
        '''
            初始化时传入RASPA输出文件的字符串
        '''
        self.output_string = output_string
        self.components = re.findall(
            r'Component \d+ \[(.*)\] \(Adsorbate molecule\)', self.output_string)

    def get_components(self):
        return self.components

    def is_finished(self):
        '''
            返回该任务是否已完成
        '''
        pattern = r'Simulation finished'
        result = re.findall(pattern, self.output_string)
        return len(result) > 0

    def get_warnings(self):
        '''
            返回存储警告信息的列表
        '''
        if len(re.findall(r'0 warnings', self.output_string)) > 0:
            return []
        pattern = r'WARNING: (.*)\n'
        return list(set(re.findall(pattern, self.output_string)))

    def get_pressure(self):
        '''
            返回压力，单位是Pa
        '''
        pattern = r'Pressure:\s+(.*)\s+\[Pa\]'
        result = re.findall(pattern, self.output_string)
        return result[0]

    def get_excess_adsorption(self, unit='cm^3/g'):
        '''
            指定单位，返回超额吸附量，返回值是一个字典，键是吸附质的名称，值是吸附量
            若不指定单位，默认为cm^3/g
            unit: 'mol/uc','cm^3/g','mol/kg','mg/g','cm^3/cm^3'
        '''
        patterns = {'mol/uc': r"Average loading excess \[molecules/unit cell\]\s+(-?\d+\.?\d*)\s+",
                    'cm^3/g': r"Average loading excess \[cm\^3 \(STP\)/gr framework\]\s+(-?\d+\.?\d*)\s+",
                    'mol/kg': r"Average loading excess \[mol/kg framework\]\s+(-?\d+\.?\d*)\s+",
                    'mg/g': r"Average loading excess \[milligram/gram framework\]\s+(-?\d+\.?\d*)\s+",
                    'cm^3/cm^3': r"Average loading excess \[cm\^3 \(STP\)/cm\^3 framework\]\s+(-?\d+\.?\d*)\s+"
                    }
        if unit not in patterns.keys():
            raise ValueError('单位错误！')
        result = {}
        data = re.findall(patterns[unit], self.output_string)
        for i, j in zip(self.components, data):
            result[i] = j
        return result

    def get_absolute_adsorption(self, unit='cm^3/g'):
        '''
            指定单位，返回绝对吸附量，返回值是一个字典，键是吸附质的名称，值是吸附量;
            若不指定单位，默认为cm^3/g
            unit: 'mol/uc','cm^3/g','mol/kg','mg/g','cm^3/cm^3'
        '''
        patterns = {'mol/uc': r"Average loading absolute \[molecules/unit cell\]\s+(-?\d+\.?\d*)\s+",
                    'cm^3/g': r"Average loading absolute \[cm\^3 \(STP\)/gr framework\]\s+(-?\d+\.?\d*)\s+",
                    'mol/kg': r"Average loading absolute \[mol/kg framework\]\s+(-?\d+\.?\d*)\s+",
                    'mg/g': r"Average loading absolute \[milligram/gram framework\]\s+(-?\d+\.?\d*)\s+",
                    'cm^3/cm^3': r"Average loading absolute \[cm\^3 \(STP\)/cm\^3 framework\]\s+(-?\d+\.?\d*)\s+"
                    }
        if unit not in patterns.keys():
            raise ValueError('单位错误！')
        result = {}
        data = re.findall(patterns[unit], self.output_string)
        for i, j in zip(self.components, data):
            result[i] = j
        return result

def get_unit_cell(cif_location, cutoff):
    with open(cif_location, 'r') as f:
        text = f.readlines()
    for i in text:
        if (i.startswith('_cell_length_a')):
            a = float(i.split()[-1].strip().split('(')[0])
        elif (i.startswith('_cell_length_b')):
            b = float(i.split()[-1].strip().split('(')[0])
        elif (i.startswith('_cell_length_c')):
            c = float(i.split()[-1].strip().split('(')[0])
        elif (i.startswith('_cell_angle_alpha')):
            alpha = float(i.split()[-1].strip().split('(')[0]) * math.pi / 180
        elif (i.startswith('_cell_angle_beta')):
            beta = float(i.split()[-1].strip().split('(')[0]) * math.pi / 180
        elif (i.startswith('_cell_angle_gamma')):
            gamma = float(i.split()[-1].strip().split('(')[0]) * math.pi / 180
            break
    #计算晶胞体积;pi = 3.1415926
    V =  a * b * c * (1 + 2 * math.cos(alpha) * math.cos(beta) * math.cos(gamma) - (math.cos(alpha))**2 - (math.cos(beta))**2 - (math.cos(gamma))**2) ** 0.5

    # 计算晶胞各个面的表面积
    base_area_x = b * c * math.sin(alpha)
    base_area_y = a * c * math.sin(beta)
    base_area_z = a * b * math.sin(gamma)

    # 计算各个方向的最小距离，即平行六面体各个方向的高，等于体积除以底面积
    perpendicular_length_x = V / base_area_x
    perpendicular_length_y = V / base_area_y
    perpendicular_length_z = V / base_area_z

    # 根据截断半径cutoff计算所需各方向的unit_cell数目
    a_unitcell = math.ceil(2 * cutoff / perpendicular_length_x)
    b_unitcell = math.ceil(2 * cutoff / perpendicular_length_y)
    c_unitcell = math.ceil(2 * cutoff / perpendicular_length_z)

    return "{} {} {}".format(a_unitcell, b_unitcell, c_unitcell)


def generate_simulation_input(template: str, temperature: str, pressure: str, cutoff: float, cif_dir: str,
                              cif_file: str):
    unitcell = get_unit_cell(os.path.join(cif_dir, cif_file), cutoff)
    cif_name = cif_file[:-4]
    return template.format(cif_name=cif_name, cutoff=cutoff, unitcell=unitcell, temperature=temperature,
                           pressure=pressure)


def work(cif_dir: str, cif_file: str, RASPA_dir: str, pressure: str, input_text: str, lock: Lock, q: Queue):
    cif_name = cif_file[:-4]
    curr_dir = os.path.abspath(os.path.dirname(__file__))
    result_file = os.path.join(curr_dir, "results", cif_name + "_result.csv")
    components = get_components_from_input(input_text)
    headers = get_field_headers(components)
    cmd_dir = os.path.join(curr_dir, "RASPA_Output", cif_name, pressure)
    lock.acquire()
    if not os.path.exists(result_file):
        with open(result_file, 'w') as f:
            for i in range(len(headers)):
                if i != len(headers) - 1:
                    f.write(headers[i] + ",")
                else:
                    f.write(headers[i] + "\n")

    if not os.path.exists(cmd_dir):
        os.makedirs(cmd_dir)
    lock.release()

    shutil.copy(os.path.join(cif_dir, cif_file), cmd_dir)
    cmd = os.path.join(RASPA_dir, "bin", "simulate") + " simulation.input"
    with open(os.path.join(cmd_dir, "simulation.input"), "w") as f1:
        f1.write(input_text)
        f1.close()
    os.chdir(cmd_dir)
    if os.system(cmd) == 0:
        lock.acquire()
        try:
            output_file = os.listdir(os.path.join(
                cmd_dir, "Output", "System_0"))[0]
            with open(os.path.join(cmd_dir, "Output", "System_0", output_file), 'r') as f2:
                result = get_result(f2.read(), components)
                f2.close()
            write_result(result_file, result, headers)
            print("\033[0;30;42m\n{}__{} has completed\n\033[0m".format(
                cif_name, pressure))
        except Exception as e:
            write_error(result_file, pressure)
            print("\033[0;37;41m\n{}__{} error: {} !\n\033[0m".format(
                cif_name, pressure, repr(e)))
        lock.release()
    else:
        lock.acquire()
        write_error(result_file, pressure)
        print("\033[0;37;41m\n{}__{} error !\n\033[0m".format(
            cif_name, pressure))
        lock.release()
    q.put(1)


def get_result(output_str: str, components: list):
    res = {}
    units = ['mol/uc', 'cm^3/g', 'mol/kg', 'mg/g', 'cm^3/cm^3']
    output = RASPA_Output_Data(output_str)
    res["pressure"] = output.get_pressure()
    res["finished"] = str(output.is_finished())
    res["warning"] = ""
    if res["finished"] == 'True':
        for w in output.get_warnings():
            res["warning"] += (w + "; ")

        for unit in units:
            absolute_capacity = output.get_absolute_adsorption(unit=unit)
            excess_capacity = output.get_excess_adsorption(unit=unit)
            for c in components:
                res[c + "_absolute_" + unit] = absolute_capacity[c]
                res[c + "_excess_" + unit] = excess_capacity[c]
    else:
        for unit in units:
            for c in components:
                res[c + "_absolute_" + unit] = " "
                res[c + "_excess_" + unit] = " "
    return res


def get_field_headers(components: list):
    headers = ["pressure", "finished"]
    units = ['mol/uc', 'cm^3/g', 'mol/kg', 'mg/g', 'cm^3/cm^3']
    for i in ["absolute", "excess"]:
        for j in components:
            for unit in units:
                headers.append(j + "_" + i + "_" + unit)
    headers.append("warning")
    return headers


def get_components_from_input(input_text: str):
    components = re.findall(r'MoleculeName\s+(.+)', input_text)
    return components


def write_result(result_file, result: dict, headers: list):
    with open(result_file, 'a') as f:
        for i in range(len(headers)):
            if i != len(headers) - 1:
                f.write(result[headers[i]] + ",")
            else:
                f.write(result[headers[i]] + "\n")
        f.close()


def write_error(result_file, pressure):
    with open(result_file, 'a') as f:
        f.write(pressure + ",Error,\n")
        f.close()


def check_parameters():
    cur_path = os.path.abspath(os.path.dirname(__file__))
    os.chdir(cur_path)
    config = configparser.ConfigParser()
    config.read("config.ini", encoding='utf8')
    section = "ISOTHERM_CONFIG"
    full_options = ['raspa_dir', 'cif_location',
                    'temperature', 'pressures', 'cutoffvdm', 'max_threads']
    options_in_config = config.options(section)
    missing_options = []
    option_dic = {}
    for op in full_options:
        if op not in options_in_config:
            missing_options.append(op)
        else:
            option_dic[op] = config.get(section, op)

    if len(missing_options) > 0:
        print("配置文件中参数不完整! (The parameters in the configuration file are incomplete !)")
        print("缺少的选项 (missing options) : " + str(missing_options))
        exit()

    raspa_dir = option_dic['raspa_dir']
    cif_dir = option_dic['cif_location']
    temperature = option_dic['temperature']
    pressures_str = option_dic['pressures']
    cutoffvdm = option_dic['cutoffvdm']
    max_threads = option_dic['max_threads']

    if len(raspa_dir) > 0:
        raspa_dir = os.path.abspath(raspa_dir)

    if len(cif_dir) > 0:
        cif_dir = os.path.abspath(cif_dir)

    if not os.path.exists(os.path.join(raspa_dir, "bin", "simulate")):
        print('RASPA目录无效！(Invalid RASPA_dir!)')
        exit()

    if not os.path.exists(cif_dir):
        print('cif目录无效！(Invalid cif_location!)')
        exit()

    pressures = pressures_str.split(",")
    try:
        for p in pressures:
            float(p)
    except:
        print("压力必须为数字或者使用科学计数法！(Pressures must be numerical or use scientific notation !)")
        exit()

    try:
        float(temperature)
    except:
        print("温度必须为数字！(Temperature must be numerical !)")
        exit()

    try:
        cutoffvdm = float(cutoffvdm)
    except:
        print("截断半径必须为数字！(CutOffVDM must be numerical !)")
        exit()

    try:
        max_threads = int(max_threads)
    except:
        print("线程数必须为整数！(max_threads must be integer !)")
        exit()

    if os.path.isfile(cif_dir):
        cifs = []
        cifs.append(os.path.basename(cif_dir))
        cif_dir = os.path.dirname(cif_dir)
        return raspa_dir, cif_dir, cifs, temperature, pressures, cutoffvdm, max_threads

    cifs = os.listdir(cif_dir)
    dels = []
    for cif in cifs:
        if not cif.endswith('.cif'):
            dels.append(cif)
    for s in dels:
        cifs.remove(s)
    if len(cifs) == 0:
        print('cif目录中缺乏有效的cif文件！(There are no valid cif files in the cif_location)')
        exit()

    return raspa_dir, cif_dir, cifs, temperature, pressures, cutoffvdm, max_threads


def main():
    cur_path = os.path.abspath(os.path.dirname(__file__))
    os.chdir(cur_path)
    raspa_dir, cif_dir, cifs, temperature, pressures, cutoffvdm, max_threads = check_parameters()
    
    # 设置环境变量(如果不设置，slurm系统可能出现raspa路径错误)
    os.environ['RASPA_DIR'] = raspa_dir
    os.environ['LD_LIBRARY_PATH'] = os.path.join(raspa_dir, "lib")
    
    lock = Lock()

    with open("./simulation_template.input", "r") as f:
        template = f.read()

    results_dir = os.path.join(cur_path, "results")
    if os.path.exists(results_dir):
        print("results目录已存在，请手动删除后重试！(The results fold already exists, please delete it and try again !)")
        exit()
    raspa_output_dir = os.path.join(cur_path, "RASPA_Output")
    if os.path.exists(raspa_output_dir):
        print("RASPA_Output 目录已存在，请手动删除后重试！(The RASPA_Output fold already exists, please delete it and try again !)")
        exit()
    os.makedirs(results_dir)
    os.makedirs(raspa_output_dir)

    q = Queue(maxsize=max_threads)
    for i in range(max_threads):
        q.put(1)
    for cif in cifs:
        for pressure in pressures:
            q.get()
            input_text = generate_simulation_input(template=template, temperature=temperature, pressure=pressure,
                                                   cutoff=cutoffvdm, cif_dir=cif_dir, cif_file=cif)
            thread = threading.Thread(target=work, args=(cif_dir, cif, raspa_dir,
                                                         pressure, input_text, lock, q))

            thread.start()
            time.sleep(0.4)
            os.chdir(cur_path)

    for t in threading.enumerate():
        if t.is_alive() and t.name != "MainThread":
            t.join()

    print("\033[0;30;42m\n完成！(Finish)\n\033[0m")


if __name__ == '__main__':
    main()
