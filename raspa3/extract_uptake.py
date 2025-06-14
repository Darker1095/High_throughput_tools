import re
import os
import csv
import glob
from collections import OrderedDict

def extract_adsorption_data(file_path):
    with open(file_path, 'r', encoding='utf-8') as file:
        content = file.read()

    # 定位"Loadings"部分
    loadings_start = content.find("Loadings")
    if loadings_start == -1:
        return None
    loadings_section = content[loadings_start:]

    # 匹配每个组分的完整部分
    comp_pattern = re.compile(
        r'Component\s+\d+\s+\((.*?)\)\s*([\s\S]*?)(?=Component\s+\d+\s+\(|\Z)',
        re.DOTALL
    )

    result = OrderedDict()
    comp_blocks = comp_pattern.findall(loadings_section)

    if not comp_blocks:
        return None

    # 匹配Abs和Excess的平均值
    values_pattern = re.compile(
        r'Abs\. loading average\s+([\d\.eE+-]+).*?\[molecules/cell\].*?'
        r'Abs\. loading average\s+([\d\.eE+-]+).*?\[mol/kg-framework\].*?'
        r'Abs\. loading average\s+([\d\.eE+-]+).*?\[mg/g-framework\][\s\S]*?'
        r'Excess loading average\s+([\d\.eE+-]+).*?\[molecules/cell\].*?'
        r'Excess loading average\s+([\d\.eE+-]+).*?\[mol/kg-framework\].*?'
        r'Excess loading average\s+([\d\.eE+-]+).*?\[mg/g-framework\]',
        re.DOTALL
    )

    for comp_name, comp_content in comp_blocks:
        match = values_pattern.search(comp_content)
        if not match:
            continue

        abs_vals = match.groups()[0:3]
        exc_vals = match.groups()[3:6]

        comp_prefix = comp_name.strip()
        result[f"{comp_prefix}_absolute_molecules/uc"] = float(abs_vals[0])
        result[f"{comp_prefix}_absolute_mol/kg"] = float(abs_vals[1])
        result[f"{comp_prefix}_absolute_mg/g"] = float(abs_vals[2])
        result[f"{comp_prefix}_excess_molecules/uc"] = float(exc_vals[0])
        result[f"{comp_prefix}_excess_mol/kg"] = float(exc_vals[1])
        result[f"{comp_prefix}_excess_mg/g"] = float(exc_vals[2])

    return result

def get_all_components_and_units(rows):
    components = set()
    units = ['absolute_molecules/uc', 'absolute_mol/kg', 'absolute_mg/g',
             'excess_molecules/uc', 'excess_mol/kg', 'excess_mg/g']
    for row in rows:
        for key in row.keys():
            if key == 'filename':
                continue
            comp = key.split('_', 1)[0]
            components.add(comp)
    components = sorted(components)
    return components, units

def write_to_csv(rows, output_file):
    if not rows:
        print("No data to write.")
        return

    components, units = get_all_components_and_units(rows)
    headers = ['name']
    for comp in components:
        for unit in units:
            headers.append(f"{comp}_{unit}")
    
    # 获取当前路径下的cif文件名
    cif_files = [f for f in os.listdir('.') if f.endswith('.cif')]
    cif_name = os.path.splitext(cif_files[0])[0]
    for row in rows:
        row['name'] = cif_name

    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()
        writer.writerows(rows)
                               
if __name__ == "__main__":
    files = glob.glob("output/output*.txt")
    results = []
    file_count, success_count = 0, 0
    for file in files:
        data = extract_adsorption_data(file)
        results.append(data)
        # 提取温度和压力
        match = re.search(r'output_(\d+)_([\deE\+\-]+)\.s\d+\.txt', os.path.basename(file))
        if match:
            temp = match.group(1)
            pressure = float(match.group(2))/1e5
            # 提取组分名
            if data:
                comps = [k.split('_')[0] for k in data.keys() if '_absolute_molecules/uc' in k]
                comps_str = '_'.join(comps)
            else:
                comps_str = 'unknown'
            csv_name = f"{comps_str}_{temp}K_{pressure:.1f}bar.csv"
        else:
            csv_name = "adsorption_loadings_summary.csv"
        write_to_csv(results, csv_name)
