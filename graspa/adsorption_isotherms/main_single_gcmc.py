import configparser
import math
import os
import re
import shutil
import threading
import time
import subprocess
from queue import Queue
from threading import Lock

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


def generate_simulation_input(template: str, cutoff: float, cif_dir: str,
                              cif_file: str, temperature: float, pressure: float):
    unitcell = get_unit_cell(os.path.join(cif_dir, cif_file), cutoff)
    cif_name = cif_file[:-4]
    return template.format(cif_name=cif_name, cutoff=cutoff, unitcell=unitcell, Temperature=temperature, Pressure=pressure)


def process_forcefield_files(cmd_dir, input_text):
    # 从输入文本中获取框架名称
    def get_frameworks_from_input(input_text):
        frameworks = []
        for line in input_text.splitlines():
            if "FrameworkName" in line:
                spline = line.split()
                spline.pop(0)
                frameworks = spline
        return frameworks

    # 从输入文本中获取吸附物名称
    def get_adsorbates_from_input(input_text):
        adsorbates = []
        for line in input_text.splitlines():
            if "Component" in line:
                spline = line.split()
                adsorbates.append(spline[3])
        return adsorbates

    # 读取框架伪原子的标签
    def read_Framework_pseudoAtoms(filename):
        ncols = 0
        count = 0
        labels = []
        with open(filename) as f:
            for line in f:
                if("_atom_" in line):
                    ncols += 1
                if(ncols > 0 and count > 0 and not("_atom_" in line)):
                    spline = line.split()
                    if not (ncols == len(spline)):
                        raise Exception("CIF FILE WRONG!")
                    label = spline[0]
                    label = ''.join(i for i in label if not i.isdigit())
                    labels.append(label)
                count += 1
        labels = list(set(labels))
        return labels

    # 读取分子伪原子的标签
    def read_Molecule_pseudoAtoms(filename):
        count = 0
        NAtom = 0
        Atomcount = 0
        labels = []
        with open(filename) as f:
            for line in f:
                if(count == 5):
                    spline = line.split()
                    NAtom = int(spline[0])
                if("atomic positions" in line):
                    start = count
                if(NAtom > 0 and 'atomic positions' in line):
                    Atomcount = 0
                if(NAtom > 0 and Atomcount < NAtom and not("#" in line) and count > 5):
                    spline = line.split()
                    if len(spline) > 1:
                        labels.append(spline[1])
                        Atomcount += 1
                count += 1
        labels = list(set(labels))
        return labels

    # 处理力场文件
    def Process_ForceFieldFile(Labels, dir_path):
        input_file    = os.path.join(dir_path, "force_field_mixing_rules.def")
        original_file = os.path.join(dir_path, "original_force_field_mixing_rules.def")
        output_file   = os.path.join(dir_path, "output_force_field_mixing_rules.def")
        shutil.copy(input_file, original_file)
        count = 0
        newlines=len(Labels)
        oldlines=0
        with open(original_file) as fr:
            with open(output_file, 'w') as fw:
                for line in fr:
                    if(count == 5):
                        spline = line.split()
                        oldlines = int(spline[0])
                        line = str(int(newlines)) + '\n'
                    if(count < 7 or count >= (7 + oldlines)):
                        fw.write(line)
                    else:
                        spline  = line.split()
                        if spline:
                            element = spline[0]
                            if(element in Labels):
                                fw.write(line)
                    count += 1

    # 处理伪原子文件
    def Process_PseudoAtomsFile(Labels, dir_path):
        input_file    = os.path.join(dir_path, "pseudo_atoms.def")
        original_file = os.path.join(dir_path, "original_pseudo_atoms.def")
        output_file   = os.path.join(dir_path, "output_pseudo_atoms.def")
        shutil.copy(input_file, original_file)
        count = 0
        newlines=len(Labels)
        oldlines=0
        with open(original_file) as fr:
            with open(output_file, 'w') as fw:
                for line in fr:
                    if(count == 1):
                        spline = line.split()
                        oldlines = int(spline[0])
                        line = str(int(newlines)) + '\n'
                    if(count < 3):
                        fw.write(line)
                    else:
                        spline  = line.split()
                        if spline:
                            element = spline[0]
                            if(element in Labels):
                                fw.write(line)
                    count += 1

    # 获取框架名称和吸附物名称
    frameworks = get_frameworks_from_input(input_text)
    adsorbates = get_adsorbates_from_input(input_text)
    Frameworklabels = []
    Adsorbatelabels = []
    # 遍历框架名称，读取框架伪原子的标签
    for framework in frameworks:
        cif_path = os.path.join(cmd_dir, framework + '.cif')
        if os.path.exists(cif_path):
            Frameworklabels.extend(read_Framework_pseudoAtoms(cif_path))
    # 遍历吸附物名称，读取分子伪原子的标签
    for adsorbate in adsorbates:
        def_path = os.path.join(cmd_dir, adsorbate + '.def')
        if os.path.exists(def_path):
            Adsorbatelabels.extend(read_Molecule_pseudoAtoms(def_path))
    # 合并框架和吸附物的标签
    Labels = []
    Labels.extend(Frameworklabels)
    Labels.extend(Adsorbatelabels)
    # 处理力场文件和伪原子文件
    Process_ForceFieldFile(Labels, cmd_dir)
    Process_PseudoAtomsFile(Labels, cmd_dir)
    # 覆盖原始文件
    shutil.move(os.path.join(cmd_dir, "output_force_field_mixing_rules.def"), os.path.join(cmd_dir, "force_field_mixing_rules.def"))
    shutil.move(os.path.join(cmd_dir, "output_pseudo_atoms.def"), os.path.join(cmd_dir, "pseudo_atoms.def"))
    # 删除 original 文件
    os.remove(os.path.join(cmd_dir, "original_force_field_mixing_rules.def"))
    os.remove(os.path.join(cmd_dir, "original_pseudo_atoms.def"))


def wait_for_task_finish(output_txt_path, timeout=3600, interval=5):
    """等待Output.txt出现END OF PROGRAM，超时单位秒"""
    waited = 0
    while waited < timeout:
        if os.path.exists(output_txt_path):
            with open(output_txt_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
                if "END OF PROGRAM" in content:
                    return True
        time.sleep(interval)
        waited += interval
    return False

def work(cif_dir: str, cif_file: str, gRASPA_dir: str, result_file: str, components: str, headers: str, input_text: str,
         lock: Lock, q: Queue, pressure: float):
    cif_name = cif_file[:-4]
    curr_dir = os.path.abspath(os.path.dirname(__file__))
    output_dir = os.path.join(curr_dir, "gRASPA_Output")
    # 新增压力子目录
    cmd_dir = os.path.join(output_dir, cif_name, str(pressure))
    if not os.path.exists(cmd_dir):
        os.makedirs(cmd_dir)
    # 复制cif文件到工作目录
    shutil.copy(os.path.join(cif_dir, cif_file), cmd_dir)
    
    # 复制single_FF目录下的所有文件到工作目录
    FF_dir = os.path.join(curr_dir, "single_FF")   
    for filename in os.listdir(FF_dir):
        src_file = os.path.join(FF_dir, filename)
        dst_file = os.path.join(cmd_dir, filename)
        if os.path.isfile(src_file):
            shutil.copy(src_file, dst_file)

    # 处理力场文件
    process_forcefield_files(cmd_dir, input_text)
    
    cmd = ["gRASPA"]
    sim_input_path = os.path.join(cmd_dir, "simulation.input")
    with open(sim_input_path, "w") as f1:
        f1.write(input_text)

    # Run gRASPA in the specific command directory without changing global cwd
    output_txt_path = os.path.join(cmd_dir, "output.txt")
    try:
        with open(output_txt_path, "w") as out:
            subprocess.run(cmd, stdout=out, stderr=subprocess.STDOUT, cwd=cmd_dir, check=False)
    except Exception:
        # If subprocess.run fails for any reason, proceed to waiting/checking as before
        pass

    if wait_for_task_finish(output_txt_path):
        lock.acquire()
        try:
            output_dir = os.path.join(cmd_dir, "Output")
            output_file = next(f for f in os.listdir(output_dir) if f.startswith("System_0") and f.endswith(".data"))
            with open(os.path.join(output_dir, output_file), 'r') as f2:
                result = get_result(f2.read(), components, pressure)
            write_result(result_file, result, headers)
            print("\033[0;30;42m\n{} has completed\n\033[0m".format(cif_name))
        except Exception as e:
            write_error(result_file, pressure)
            print("\033[0;37;41m\n{} error: {} !\n\033[0m".format(cif_name, repr(e)))
        lock.release()
    else:
        lock.acquire()
        write_error(result_file, pressure)
        print("\033[0;37;41m\n{} error: Output.txt not finished!\n\033[0m".format(cif_name))
        lock.release()
    q.put(1)


def get_result(output_str: str, components: list, pressure: float):
    # 单组分的吸附数据
    content = output_str.splitlines()
    res = {}
    res["pressure"] = str(pressure)
    res["finished"] = "True"
    res["warning"] = ""

    avg_heat_of_adsorption = ""
    loading_units = {
        "# MOLECULES": "",
        'Framework mass': "",
        "mg/g": "",
        "mol/kg": "",
        "g/L": ""
    }

    for i, line in enumerate(content):
        if "BLOCK AVERAGES (HEAT OF ADSORPTION: kJ/mol)" in line:
            try:
                avg_heat_of_adsorption = content[i+7].split(',')[0].split()[-1].strip()
            except Exception:
                pass
        if 'BLOCK AVERAGES (LOADING: # MOLECULES)' in line:
            try:
                loading_units["# MOLECULES"] = content[i+16].split(',')[0].split()[-1].strip()
            except Exception:
                pass
        if 'BLOCK AVERAGES (LOADING: mg/g)' in line:
            try:
                loading_units['Framework mass'] = content[i+3].split(' ')[3].strip()
                loading_units["mg/g"] = content[i+19].split(',')[0].split()[-1].strip()
            except Exception:
                pass
        if 'BLOCK AVERAGES (LOADING: mol/kg)' in line:
            try:
                loading_units["mol/kg"] = content[i+19].split(',')[0].split()[-1].strip()
            except Exception:
                pass
        if 'BLOCK AVERAGES (LOADING: g/L)' in line:
            try:
                loading_units["g/L"] = content[i+8].split(',')[0].split()[-1].strip()
            except Exception:
                pass

    res[components[0] + "_Heat_of_adsorption_kJ/mol"] = avg_heat_of_adsorption
    res[components[0] + "_loading_molecules"] = loading_units["# MOLECULES"]
    res[components[0] + "_loading_mg/g"] = loading_units["mg/g"]
    res[components[0] + "_loading_mol/kg"] = loading_units["mol/kg"]
    res[components[0] + "_loading_g/L"] = loading_units["g/L"]

    return res


def get_field_headers(components: list):
    headers = ["pressure", "finished"]
    for c in components:
        headers.append(c + "_Heat_of_adsorption_kJ/mol")
        headers.append(c + "_loading_molecules")
        headers.append(c + "_loading_mg/g")
        headers.append(c + "_loading_mol/kg")
        headers.append(c + "_loading_g/L")          
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
        f.write(str(pressure) + ",Error,\n")
        f.close()


def check_parameters():
    cur_path = os.path.abspath(os.path.dirname(__file__))
    os.chdir(cur_path)
    config = configparser.ConfigParser()
    config.read("config.ini", encoding='utf8') 
    section = "ADSORPTION_CONFIG"
    full_options = ['gRASPA_dir', 'cif_location', 'CutOffVDM', 'max_tasks', 'Temperature', 'Pressure']
    options_in_config = [option.lower() for option in config.options(section)]
    missing_options = []
    option_dic = {}
    for op in full_options:
        if op.lower() not in options_in_config:
            missing_options.append(op)
        else:
            # 获取原始大小写的选项名
            original_case_option = next(option for option in config.options(section) if option.lower() == op.lower())
            option_dic[op] = config.get(section, original_case_option)

    if len(missing_options) > 0:
        print("配置文件中参数不完整! (The parameters in the configuration file are incomplete !)")
        print("缺少的选项 (missing options) : " + str(missing_options))
        exit()

    graspa_dir = option_dic['gRASPA_dir']
    cif_dir = option_dic['cif_location']
    cutoffvdm = option_dic['CutOffVDM']
    max_tasks = option_dic['max_tasks']
    temperature = option_dic['Temperature']
    pressure = option_dic['Pressure']

    if len(graspa_dir) > 0:
        graspa_dir = os.path.abspath(graspa_dir)

    if len(cif_dir) > 0:
        cif_dir = os.path.abspath(cif_dir)

    if not os.path.exists(graspa_dir):
        print('gRASPA目录无效！(Invalid gRASPA_dir!)')
        exit()

    if not os.path.exists(cif_dir):
        print('cif目录无效！(Invalid cif_location!)')
        exit()

    try:
        cutoffvdm = float(cutoffvdm)
    except:
        print("截断半径必须为数字！(CutOffVDM must be numerical !)")
        exit()

    try:
        max_tasks = int(max_tasks)
    except:
        print("任务数必须为整数！(max_tasks must be integer !)")
        exit()
    
    try:
        temperature = float(temperature)
    except:
        print("温度必须为数字！(Temperature must be numerical !)")
        exit()
    
    try:
        # 支持多个压力
        pressure = [float(p) for p in str(pressure).split(',')]
    except:
        print("压力必须为数字或逗号分隔的数字！(Pressure must be numerical or comma separated!)")
        exit()

    if os.path.isfile(cif_dir):
        cifs = []
        cifs.append(os.path.basename(cif_dir))
        cif_dir = os.path.dirname(cif_dir)
        return graspa_dir, cif_dir, cifs, cutoffvdm, max_tasks

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

    return graspa_dir, cif_dir, cifs, cutoffvdm, max_tasks, temperature, pressure


def main():
    cur_path = os.path.abspath(os.path.dirname(__file__))
    os.chdir(cur_path)
    graspa_dir, cif_dir, cifs, cutoffvdm, max_tasks, temperature, pressures = check_parameters()
    lock = Lock()

    with open("simulation_template.input", "r") as f:
        template = f.read()                                         

    components = get_components_from_input(template)
    if len(components) != 1:
        print("只支持单组分吸附！(Only single component adsorption is supported!)")
        exit()
    headers = get_field_headers(components)

    output_dir = os.path.join(cur_path, "gRASPA_Output")
    if os.path.exists(output_dir):
        print("gRASPA_Output目录已存在，请手动删除后重试！(The gRASPA_Output fold already exists, please delete it and try again !)")
        exit()
    os.makedirs(output_dir)

    q = Queue(maxsize=max_tasks)
    for i in range(max_tasks):
        q.put(1)

    threads = []
    for cif in cifs:
        cif_name = cif[:-4]
        # 每个cif一个csv
        result_file = os.path.join(cur_path, f"{cif_name}.csv")
        if os.path.exists(result_file):
            os.remove(result_file)
        with open(result_file, 'w') as f:
            for i in range(len(headers)):
                if i != len(headers) - 1:
                    f.write(headers[i] + ",")
                else:
                    f.write(headers[i] + "\n")
            f.close()
        for p in pressures:
            q.get()
            input_text = generate_simulation_input(
                template=template, cutoff=cutoffvdm, cif_dir=cif_dir, cif_file=cif, 
                temperature=float(temperature), pressure=float(p))
            thread = threading.Thread(target=work, args=(cif_dir, cif, graspa_dir,
                                                         result_file, components, headers, input_text, lock, q, p))
            thread.start()
            threads.append(thread)
            time.sleep(0.3)
            os.chdir(cur_path)

    for t in threads:
        t.join()

    print("\033[0;30;42m\n完成！(Finish)\n\033[0m".encode("utf-8").decode("latin1"))


if __name__ == '__main__':
    main()
