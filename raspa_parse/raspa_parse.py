import re

'''
示例：
    with open('./output.data','r') as f:
        str = f.read()
    parser = RASPA_Output_Data(str)
    print(parser.is_finished())
    print(parser.get_absolute_adsorption())

'''


class RASPA_Output_Data():
    '''
        RASPA输出文件对象
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
        
    def get_He_void_fraction(self):
        '''
        返回[helium] Average Widom Rosenbluth-weight字符后的数值
        '''
        pattern = r'Average Widom Rosenbluth-weight:\s+(-?\d+\.?\d*)\s+'
        result = re.findall(pattern, self.output_string)
        return result

    def get_Surface_Area(self,unit='m^2/cm^3'):
        '''
        返回Average surface area字符后的数值
        ''' 
        patterns = {'A^2':     r'Average surface area:\s+(-?\d+\.?\d*)\s+\+/-\s+(-?\d+\.?\d*)\s+\[A\^2\]', 
                   'm^2/g':                          '\s+(-?\d+\.?\d*)\s+\+/-\s+(-?\d+\.?\d*)\s+\[m\^2/g\]',
                   'm^2/cm^3':                       '\s+(-?\d+\.?\d*)\s+\+/-\s+(-?\d+\.?\d*)\s+\[m\^2/cm\^3\]'
                   }
        if unit not in patterns.keys():
            raise ValueError('单位错误！')
        result = {}
        data = re.findall(patterns[unit], self.output_string)
        result = data[0] if data else None
        return result   

    def get_adsorption_heat(self):
        '''
            返回吸附热(KJ/mol)，该数值使用波动法计算 fluctuation formula
            返回值是一个字典，键是吸附质的名称，值是吸附热;
        '''
        result = {}
        # 定义第一种情况下的正则表达式模式
        pattern1 = r'Enthalpy of adsorption component \d+ \[(.*)\]\n\s*-*\n.*\n.*\n.*\n.*\n.*\n\s*-*\n.*\n\s+(\-?\d+\.?\d*)\s+'

        # 定义第二种情况下的正则表达式模式
        pattern2 = r'Total enthalpy of adsorption\n.*\n.*\n.*\n.*\n.*\n.*\n.*\n.*\n\s+(\-?\d+\.?\d*)\s+'

        # 尝试匹配pattern1
        data1 = re.findall(pattern1, self.output_string)

        if data1:
            for i, j in zip(self.components, data1):
                result[i] = j[1]  # 使用元组中的第二个元素作为吸附热值
        else:
            # 如果pattern1匹配不成功，则匹配pattern2
            data2 = re.findall(pattern2, self.output_string)
            if data2:
                result["Total enthalpy of adsorption"] = data2[0]
        return result

    def get_adsorption_heat_infinite_dilution(self):
        '''
            返回无限稀释吸附热(KJ/mol)
            返回值是一个字典，键是吸附质的名称，值是吸附热;
        '''
        pattern = r'Total energy:\n.*\n.*\n.*\n.*\n.*\n.*\n.*\n\s+Average\s+(\-?\d+\.?\d*)\s+'
        data = re.findall(pattern, self.output_string)
        '''∆H = ∆U − RT = [Uhg] − [Uh] − [Ug] − RT
           利用该公式进行吸附热换算，∆H单位为K，框架为刚性Uh = 0，气体分子能量Ug=0,能量主要来自气体分子与框架的相互作用
        '''
        data = [float(x) for x in data]
        result = (data[0] - 300) * 8.314462618/1000
        return result

    def get_henry_coefficient(self):
        '''
            返回亨利系数(mol/kg/Pa)
            返回值是一个字典，键是吸附质的名称，值是亨利系数;
        '''
        pattern = r'\[.*\]\s+Average Henry coefficient:\s+(-?\d+\.\d+e[+-]\d+|-?\d+\.?\d*)\s+'
        data = re.findall(pattern, self.output_string)
        result = {}
        for i, j in zip(self.components, data):
            result[i] = j
        return result
    
    def get_Framework_density(self):
        '''
        返回Framework Density字符后的数值，单位是kg/m^3
        '''
        pattern = r'Framework Density:\s+(-?\d+\.?\d*)\s+\[kg/m\^3\]\s+'
        result = re.findall(pattern, self.output_string)
        return result
    
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

    def get_adsorption_heat(self):
        '''
            返回吸附热(KJ/mol)
            返回值是一个字典，键是吸附质的名称，值是吸附热;
        '''
        result = {}
        if len(self.components) > 1:
            pattern = r'Component \d+ \[(.*)\]\n\s*-*\n.*\n.*\n.*\n.*\n.*\n\s*-*\n.*\n\s+(\-?\d+\.?\d*)\s'
        else:
            pattern = r'Enthalpy of adsorption:\n.*\n.*\n.*\n.*\n.*\n.*\n.*\n.*\n.*\n\s+(\-?\d+\.?\d*)\s'
        data = re.findall(pattern, self.output_string)
        for i, j in zip(self.components, data):
            result[i] = j
        return result

    def get_henry_coefficient(self):
        '''
            返回亨利系数(mol/kg/Pa)
            返回值是一个字典，键是吸附质的名称，值是亨利系数;
        '''
        pattern = r'\[.*\]\s+Average Henry coefficient:\s+(-?\d+\.?\d*)\s+'
        data = re.findall(pattern, self.output_string)
        result = {}
        for i, j in zip(self.components, data):
            result[i] = j
        return result

    def get_all_adsorption_result(self):
        '''
            返回值是一个字典 (dic)，包含各组分各单位的吸附量数据
            键的命名格式为"component_absulute_unit" (绝对吸附量)，"component_absulute_unit" (超额吸附量)
            例如dic["H2_absolute_mol/kg"]的值表示H2的绝对吸附量，单位是mol/kg
            此外，dic["finished"]表示是否已完成，dic["warning"]表示警告信息
        '''
        res = {}
        units = ['mol/uc', 'cm^3/g', 'mol/kg', 'mg/g', 'cm^3/cm^3']
        res["finished"] = str(self.is_finished())
        res["warning"] = ""
        if res["finished"] == 'True':
            for w in self.get_warnings():
                res["warning"] += (w + "; ")

            for unit in units:
                absolute_capacity = self.get_absolute_adsorption(unit=unit)
                excess_capacity = self.get_excess_adsorption(unit=unit)
                for c in self.components:
                    res[c + "_absolute_" + unit] = absolute_capacity[c]
                    res[c + "_excess_" + unit] = excess_capacity[c]
        else:
            for unit in units:
                for c in self.components:
                    res[c + "_absolute_" + unit] = " "
                    res[c + "_excess_" + unit] = " "
        return res
