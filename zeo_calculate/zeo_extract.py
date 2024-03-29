import os
import csv

def get_LCD_PLD(string):
    strs = string.split()
    LCD = strs[3]
    PLD = strs[2]
    return {'LCD': LCD, 'PLD': PLD}

def get_density_VSA_GSA(string):
    strs = string.split()
    density = strs[5]
    VSA = strs[9]
    GSA = strs[11]
    return {'density(g/cm^3)': density, 'VSA(m^2/cm^3)': VSA, 'GSA(m^2/g)': GSA}

def get_Vp_voidFraction(string):
    strs = string.split()
    Vp = strs[11]
    void_fraction = strs[9]
    return {'Vp(cm^3/g)': Vp, 'void_fraction': void_fraction}

def extract_info_from_file(file_path):
    with open(file_path, 'r') as file:
        content = file.read()
        if file_path.endswith('.res'):
            return get_LCD_PLD(content)
        elif file_path.endswith('.sa'):
            return get_density_VSA_GSA(content)
        elif file_path.endswith('.vol'):
            return get_Vp_voidFraction(content)
        else:
            return {}

def main():
    cif_dir = 'cif_dir = /home/luxiuyang/RASPA_tools/test_cifs'
    results_dir = 'zeo_results'
    results = {}

    for cif_file in os.listdir(cif_dir):
        if cif_file.endswith('.cif'):
            base_name = cif_file[:-4]
            for ext in ['.res', '.sa', '.vol']:
                file_path = os.path.join(results_dir, base_name + ext)
                if os.path.exists(file_path):
                    info = extract_info_from_file(file_path)
                    if base_name not in results:
                        results[base_name] = info
                    else:
                        results[base_name].update(info)

    with open('myresults.csv', 'w', newline='') as csvfile:
        fieldnames = ['name', 'LCD', 'PLD', 'density(g/cm^3)', 'VSA(m^2/cm^3)', 'GSA(m^2/g)', 'Vp(cm^3/g)', 'void_fraction']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for name, info in results.items():
            info['name'] = name
            writer.writerow(info)

if __name__ == "__main__":
    main()