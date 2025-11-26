# 添加在文件最顶部（import之前）
import sys
import os

# 解决PyInstaller路径问题
if getattr(sys, 'frozen', False):
    # 打包后的路径
    BASE_DIR = sys._MEIPASS
else:
    # 开发环境路径
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# 确保RDKit能正确找到数据文件
os.environ['RDKIT_DATA_DIR'] = os.path.join(BASE_DIR, 'rdkit/Data')

import sys
import numpy as np
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout,
                             QHBoxLayout, QLineEdit, QPushButton, QLabel,
                             QComboBox, QFrame, QMessageBox, QSizePolicy,
                             QTableWidget, QTableWidgetItem, QHeaderView,
                             QCheckBox, QRadioButton, QButtonGroup)
from PyQt5.QtCore import Qt, QTimer
from PyQt5.QtGui import QFont, QBrush, QColor, QPainter, QPen, QBrush
from PyQt5.QtOpenGL import QGLWidget
from OpenGL.GL import *
from OpenGL.GLU import *
import logging
import webbrowser
import math

# 配置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def check_numpy_version():
    """检查NumPy版本是否与RDKit兼容"""
    try:
        import numpy as np
        version = np.__version__
        # 检查是否为NumPy 2.x
        if version.startswith('2.'):
            return False, version
        return True, version
    except Exception as e:
        logger.error(f"Error checking NumPy version: {str(e)}")
        return False, "unknown"


def check_dependencies():
    """检查必要的依赖是否安装"""
    missing_deps = []

    # 检查NumPy版本
    numpy_compatible, numpy_version = check_numpy_version()
    if not numpy_compatible:
        missing_deps.append(f"NumPy版本不兼容 (当前: {numpy_version})。RDKit需要NumPy 1.x。\n"
                            "请运行: pip install 'numpy<2'")

    # 检查PyQt5
    try:
        from PyQt5.QtWidgets import QApplication
    except ImportError:
        missing_deps.append("PyQt5 (pip install pyqt5)")

    # 检查PyOpenGL
    try:
        from OpenGL.GL import glClear
    except ImportError:
        missing_deps.append("PyOpenGL (pip install PyOpenGL PyOpenGL_accelerate)")

    # 检查RDKit
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, Descriptors, Lipinski, rdMolDescriptors
        # 确保我们可以访问所有必要的函数
        _ = rdMolDescriptors.CalcMolFormula
        _ = rdMolDescriptors.CalcTPSA
    except ImportError:
        missing_deps.append("RDKit (conda install -c conda-forge rdkit 或 pip install rdkit-pypi)")
    except Exception as e:
        missing_deps.append(f"RDKit加载错误: {str(e)}\n"
                            "这通常是NumPy版本不兼容导致的。请确保安装了NumPy 1.x版本\n"
                            "运行: pip install 'numpy<2'")

    return missing_deps


# 在导入RDKit前检查NumPy版本
numpy_compatible, numpy_version = check_numpy_version()
if not numpy_compatible:
    logger.warning(f"NumPy version {numpy_version} is not compatible with RDKit. RDKit requires NumPy 1.x.")
else:
    # 只有在NumPy版本兼容时才导入RDKit
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, Descriptors, Lipinski, rdMolDescriptors

        logger.info("RDKit imported successfully")
    except Exception as e:
        logger.error(f"Failed to import RDKit even with compatible NumPy: {str(e)}")

# 常见分子名称到SMILES的映射（已修复语法错误）
MOLECULE_DB = {
    "水": "O",
    "苯": "c1ccccc1",
    "水分子": "O",
    "H2O": "O",
    "氨": "N",
    "氨气": "N",
    "NH3": "N",
    "甲烷": "C",
    "CH4": "C",
    "乙烷": "CC",
    "C2H6": "CC",
    "乙烯": "C=C",
    "C2H4": "C=C",
    "乙炔": "C#C",
    "C2H2": "C#C",
    "乙醇": "CCO",
    "酒精": "CCO",
    "C2H5OH": "CCO",
    "甲醇": "CO",
    "CH3OH": "CO",
    "丙酮": "CC(C)=O",
    "甲醛": "C=O",
    "乙醛": "CC=O",
    "丙醛": "CCC=O",
    "甲酸": "C(=O)O",
    "乙酸": "CC(=O)O",
    "醋酸": "CC(=O)O",
    "乙酸乙酯": "CC(=O)OCC",
    "丙烯酸": "C=CC(=O)O",
    "丙烯腈": "C=CC#N",
    "丙烯": "C=CC",
    "丁烷": "CCCC",
    "异丁烷": "CC(C)C",
    "丁烯": "C=CCC",
    "顺丁烯二酸": "O=C(O)/C=C/C(=O)O",
    "反丁烯二酸": "O=C(O)/C=C/C(=O)O",
    "草酸": "O=C(O)C(=O)O",
    "甘油": "C(CO)CO",
    "丙三醇": "C(CO)CO",
    "甲胺": "CN",
    "乙胺": "CCN",
    "丙胺": "CCCN",
    "甲酸甲酯": "COC(=O)",
    "甲酸乙酯": "CCOC(=O)",
    "乙酸甲酯": "COC(=O)C",
    "乙酸丙酯": "CCCOC(=O)C",
    "甲苯": "c1ccccc1C",
    "二甲苯": "c1ccc(C)cc1C",
    "乙苯": "c1ccc(CC)cc1",
    "苯乙烯": "c1ccc(C=C)cc1",
    "苯酚": "c1ccc(O)cc1",
    "甲酚": "c1ccc(C)cc1O",
    "硝基苯": "c1ccc(cc1)[N+](=O)[O-]",
    "氯苯": "c1ccc(Cl)cc1",
    "溴苯": "c1ccc(Br)cc1",
    "碘苯": "c1ccc(I)cc1",
    "苯胺": "c1ccc(N)cc1",
    "甲醚": "COC",
    "乙醚": "CCOCC",
    "丙醚": "CCCOCC",
    "硫化氢": "S",
    "H2S": "S",
    "氰化氢": "N#C",
    "氢氰酸": "N#C",
    "氯化氢": "Cl",
    "HCl": "Cl",
    "氟化氢": "F",
    "HF": "F",
    "溴化氢": "Br",
    "HBr": "Br",
    "碘化氢": "I",
    "HI": "I",
    "硫酸": "OS(=O)(=O)O",
    "硝酸": "O[N+](=O)O",
    "磷酸": "OP(=O)(O)O",
    "盐酸": "Cl",
    "醋酸酐": "CC(=O)OC(=O)C",
    "乙酸酐": "CC(=O)OC(=O)C",
    "丙酸酐": "CCC(=O)OC(=O)CC",
    "乙酰氯": "CC(=O)Cl",
    "甲酰氯": "C(=O)Cl",
    "乙酰胺": "CC(=O)N",
    "甲酰胺": "C(=O)N",
    "丙酰胺": "CCC(=O)N",
    "尿素": "C(=O)(N)N",
    "丙烯酰胺": "C=CC(=O)N",
    "丙烯酸甲酯": "C=CC(=O)OC",
    "丙烯酸乙酯": "C=CC(=O)OCC",
    "丙烯酸丁酯": "C=CC(=O)OCCCC",
    "甲基丙烯酸甲酯": "C=C(C)C(=O)OC",
    "丁二烯": "C=CCC=C",
    "异戊二烯": "C=CC(=C)C",
    "丙二酸": "C(C(=O)O)C(=O)O",
    "丁二酸": "C(CCC(=O)O)C(=O)O",
    "己二酸": "C(CCCCC(=O)O)C(=O)O",
    "对苯二甲酸": "O=C(O)c1ccc(cc1)C(=O)O",
    "间苯二甲二酸": "O=C(O)c1cccc(c1)C(=O)O",
    "邻苯二甲酸": "O=C(O)c1ccccc1C(=O)O",
    "苯甲酸": "O=C(O)c1ccccc1",
    "乙苯酸": "O=C(O)CCc1ccccc1",
    "肉桂酸": "O=C(O)/C=C/c1ccccc1",
    "乳酸": "C[C@H](O)C(=O)O",
    "苹果酸": "C([C@@H](C(=O)O)O)C(=O)O",
    "柠檬酸": "C(C(CC(=O)O)(C(=O)O)O)C(=O)O",
    "酒石酸": "C([C@@H](O)C(=O)O)(C(=O)O)O",
    "草酰乙酸": "O=C(O)CC(=O)C(=O)O",
    "α-酮戊二酸": "O=C(O)CCC(=O)C(=O)O",
    "琥珀酸": "C(CCC(=O)O)C(=O)O",
    "延胡索酸": "O=C(O)/C=C/C(=O)O",
    "马来酸": "O=C(O)/C=C/C(=O)O",
    "丙酮酸": "C(=O)C(=O)O",
    "乙醛酸": "O=CC(=O)O",
    "甘氨酸": "N[C@@H](C(=O)O)C",
    "丙氨酸": "N[C@@H](C(=O)O)CC",
    "缬氨酸": "N[C@@H](C(=O)O)C(C)C",
    "亮氨酸": "N[C@@H](C(=O)O)CC(C)C",
    "异亮氨酸": "N[C@@H](C(=O)O)[C@@H](C)CC",
    "苯丙氨酸": "N[C@@H](C(=O)O)Cc1ccccc1",
    "酪氨酸": "N[C@@H](C(=O)O)Cc1ccc(O)cc1",
    "色氨酸": "N[C@@H](C(=O)O)Cc1c[nH]c2c1cccc2",
    "丝氨酸": "N[C@@H](C(=O)O)CO",
    "苏氨酸": "N[C@@H](C(=O)O)[C@@H](O)C",
    "半胱氨酸": "N[C@@H](C(=O)O)CS",
    "蛋氨酸": "N[C@@H](C(=O)O)CCCSC",
    "脯氨酸": "N1[C@@H](CC1)C(=O)O",
    "羟脯氨酸": "N1[C@@H]([C@H](O)CC1)C(=O)O",
    "组氨酸": "N[C@@H](C(=O)O)Cc1c[nH]cn1",
    "赖氨酸": "N[C@@H](CCCNC)C(=O)O",
    "精氨酸": "N[C@@H](CCCNC(N)=N)C(=O)O",
    "天冬氨酸": "N[C@@H](CC(=O)O)C(=O)O",
    "谷氨酸": "N[C@@H](CCC(=O)O)C(=O)O",
    "天冬酰胺": "N[C@@H](CC(=O)N)C(=O)O",
    "谷氨酰胺": "N[C@@H](CCC(=O)N)C(=O)O",
    "二氧化碳": "O=C=O",
    "二氧化硫": "O=S=O",
    "三氟化硼": "FB(F)F",
    "硝酸根离子": "[N+](=O)([O-])[O-]",
    "硫酸根离子": "O=S(=O)([O-])[O-]",
    "四氟化硫": "FS(F)(F)F",
    "五氯化磷": "ClP(Cl)(Cl)(Cl)Cl",
    "六氟化硫": "S(F)(F)(F)(F)(F)F",
    "四氟化氙": "F[Xe](F)(F)F",
    "氯化铍": "Cl[Be]Cl",
    "三氧化硫": "O=S(=O)=O",
    "二溴化锡": "Br[Sn]Br",
    "二氯化铅": "Cl[Pb]Cl",
    "碳酸根": "[O-]C(=O)[O-]",
    "硝酸根": "[O-][N+](=O)=O",
    "硅酸根": "[O-][Si]([O-])([O-])[O-]",
    "四氯化碳": "C(Cl)(Cl)(Cl)Cl",
    "铵根": "[NH4+]",
    "三氯化磷": "ClP(Cl)Cl",
    "亚硫酸根": "[O-]S(=O)[O-]",
    "水合氢离子": "[OH3+]",
    "三氟化溴": "FBr(F)F",
    "二氟化氙": "FXeF",
    "亚氨根": "[NH2-]",
    "四氯化硅": "Cl[Si](Cl)(Cl)Cl",
    "碘三负离子": "[I-][I+][I-]",
    "碘三正离子": "[I+][I][I+]",
    "光气": "ClC(=O)Cl",
    "二氯亚砜": "ClS(=O)Cl",
    "氧化二氯": "ClOCl",
    "高氯酸根": "[O-]Cl(=O)(=O)=O",
    "二氧化氯": "O=[Cl]=O",
    "二氟化氧": "OFF",
    "四氟化氧": "OF(F)(F)F",
    "葡萄糖": "OC[C@H]([C@@H](O)[C@H](O)[C@H](O)CO)O",
    "氧化钙": "[Ca]=O",
    "氢氧化钠": "[Na+].[O-][H]",
    "氢氧化钙": "[Ca++](.[O-][H])(.[O-][H])",
    "氯化钠": "[Na+].[Cl-]",
    "碳酸钠": "[Na+].[Na+].C(=O)([O-])[O-]",
    "碳酸氢钠": "[Na+].C(=O)([O-])O",
    "氢气": "[H][H]",  # 单质 | 实验室制法：Zn + H₂SO₄ → H₂↑；检验：点燃淡蓝色火焰，干冷烧杯壁有水珠
    "氧气": "O=O",  # 单质 | 实验室制法：2KClO₃ → 2KCl + 3O₂↑；检验：使带火星木条复燃
    "氮气": "N#N",  # 单质 | 空气中含量 78%，稳定性源于 N≡N 三键高键能
    "氯气": "ClCl",  # 单质 | 黄绿色有毒气体；实验室制法：MnO₂ + 4HCl(浓) → Cl₂↑；检验：湿润淀粉-KI 试纸变蓝
    "臭氧": "O[O]O",  # 单质 | 淡蓝色气体，强氧化性漂白剂；与 O₂ 互为同素异形体
    "一氧化碳": "C#O",  # 氧化物 | 无色剧毒气体；还原性：Fe₂O₃ + 3CO → 2Fe + 3CO₂；与 CO₂ 区分：燃烧或通过灼热 CuO
    "一氧化氮": "[N]=O",  # 氧化物 | 无色气体，空气中氧化为红棕色 NO₂；实验室制法：3Cu + 8HNO₃(稀) → 2NO↑
    "二氧化氮": "[O][N]=O",
    "XeF₄": "FXe(F)(F)F",
    "BrF₅": "FB(F)(F)(F)F",
    "PCl₅":"ClP(Cl)(Cl)(Cl)Cl",
    # 新增：正八面体结构关键词
    "正八面体": "FS(F)(F)(F)(F)F",
    "八面体结构": "FS(F)(F)(F)(F)F",
    "octahedral": "FS(F)(F)(F)(F)F"
}

# 原始已有的颜色定义
ATOM_COLORS = {
    'H': (1.0, 1.0, 1.0),  # 白色
    'C': (0.565, 0.565, 0.565),  # 灰色
    'O': (1.0, 0.0, 0.0),  # 红色
    'N': (0.0, 0.0, 1.0),  # 蓝色
    'S': (1.0, 1.0, 0.0),  # 黄色
    'P': (1.0, 0.5, 0.0),  # 橙色
    'Cl': (0.0, 0.7, 0.0),  # 绿色
    'F': (0.0, 1.0, 1.0),  # 青色
    'Br': (0.6, 0.0, 0.0),  # 深红
    'I': (0.4, 0.0, 0.4),  # 紫色
    'Na': (0.6, 0.0, 0.8),  # 紫罗兰
    'K': (0.5, 0.3, 0.1),  # 棕色
    'Mg': (0.8, 0.8, 0.0),  # 亮黄
    'Ca': (0.6, 0.6, 0.6),  # 银灰
    'Fe': (0.8, 0.4, 0.2),  # 棕红
    'Na+': (0.6, 0.0, 0.8),  # 钠离子
    'K+': (0.5, 0.3, 0.1),  # 钾离子
    'Ca2+': (0.6, 0.6, 0.6),  # 钙离子
    'Mg2+': (0.8, 0.8, 0.0),  # 镁离子
    'Cl-': (0.0, 0.7, 0.0),  # 氯离子
    'OH-': (1.0, 0.5, 0.0)  # 氢氧根
}

# 新增的元素及其颜色
ATOM_COLORS.update({
    'Be': (0.75, 0.75, 0.75),  # 铍 - 新增，建议使用浅灰色
    'B': (1.0, 0.7, 0.7),  # 硼 - 新增，建议使用粉红色/浅红色
    'Si': (0.7, 0.7, 0.7),  # 硅 - 新增，建议使用暗灰色
    'Sn': (0.6, 0.6, 0.7),  # 锡 - 新增，建议使用蓝灰色
    'Pb': (0.5, 0.5, 0.6),  # 铅 - 新增，建议使用深灰蓝色
    'Xe': (0.3, 0.7, 0.7)  # 氙 - 新增，建议使用青绿色
})


def calculate_molecular_properties(mol):
    """计算分子的各种属性"""
    if not mol:
        return {}

    properties = {}

    try:
        # 确保分子是"干净"的 - 先进行分子规范化
        try:
            Chem.SanitizeMol(mol)
        except:
            # 如果规范化失败，尝试修复分子
            pass

        # 确保隐式价态已计算
        for atom in mol.GetAtoms():
            atom.GetImplicitValence()

        # 确保分子有氢原子
        mol_h = Chem.AddHs(mol)

        # 计算分子量
        properties["分子量"] = f"{Descriptors.ExactMolWt(mol_h):.2f} g/mol"

        # 获取分子式
        properties["分子式"] = rdMolDescriptors.CalcMolFormula(mol_h)

        # 原子和键数量
        properties["原子数量"] = mol.GetNumAtoms()
        properties["键数量"] = mol.GetNumBonds()

        # 极性判断
        polar_elements = ['O', 'N', 'F', 'Cl', 'Br', 'I']
        has_polar_bonds = False
        for bond in mol.GetBonds():
            atom1 = mol.GetAtomWithIdx(bond.GetBeginAtomIdx())
            atom2 = mol.GetAtomWithIdx(bond.GetEndAtomIdx())
            if (atom1.GetSymbol() in polar_elements and atom2.GetSymbol() not in ['C', 'H']) or \
                    (atom2.GetSymbol() in polar_elements and atom1.GetSymbol() not in ['C', 'H']):
                has_polar_bonds = True
                break

        properties["是否极性"] = "是" if has_polar_bonds else "否"

        # 计算摩尔折射率
        try:
            properties["摩尔折射率"] = f"{Chem.Descriptors.MolMR(mol_h):.2f}"
        except:
            properties["摩尔折射率"] = "计算失败"

        # 计算logP（辛醇-水分配系数）
        try:
            properties["logP"] = f"{Chem.Descriptors.MolLogP(mol_h):.2f}"
        except:
            properties["logP"] = "计算失败"

        # 计算氢键供体和受体数量
        properties["氢键供体"] = Lipinski.NumHDonors(mol_h)
        properties["氢键受体"] = Lipinski.NumHAcceptors(mol_h)

        # 计算拓扑表面积
        try:
            properties["拓扑表面积"] = f"{rdMolDescriptors.CalcTPSA(mol_h):.2f} Å²"
        except:
            properties["拓扑表面积"] = "计算失败"

        # 计算旋转键数量
        properties["旋转键数量"] = Lipinski.NumRotatableBonds(mol_h)

    except Exception as e:
        logger.error(f"Error calculating molecular properties: {str(e)}")

    return properties


class MoleculeViewer(QGLWidget):
    def __init__(self, parent=None):
        super(MoleculeViewer, self).__init__(parent)
        self.setMinimumSize(600, 600)

        # 旋转角度
        self.xRot = 0
        self.yRot = 0
        self.zRot = 0

        # 平移
        self.xTrans = 0
        self.yTrans = 0
        self.zTrans = -5  # 初始缩放

        # 鼠标位置
        self.lastX = 0
        self.lastY = 0

        # 分子数据
        self.atoms = []
        self.bonds = []
        self.atom_types = []

        # 新增：控制是否显示原子标签
        self.show_atom_labels = True
        self.show_hydrogen_labels = False

        # 标记是否已初始化
        self.initialized = False

        # 存储原子屏幕位置（用于标签渲染）
        self.atom_screen_positions = []

        # 新增：用于绘制圆柱体的显示列表
        self.bond_cylinders = {}

    def initializeGL(self):
        try:
            logger.info("Initializing OpenGL...")

            # 检查OpenGL上下文是否可用
            if not self.isValid():
                logger.error("OpenGL context is not valid")
                return

            # 设置背景颜色（更柔和的浅灰）
            glClearColor(0.85, 0.85, 0.85, 1.0)

            # 启用深度测试
            glEnable(GL_DEPTH_TEST)
            glDepthFunc(GL_LEQUAL)

            # 启用光照
            glEnable(GL_LIGHTING)
            glEnable(GL_LIGHT0)
            glEnable(GL_LIGHT1)  # 新增：第二个光源

            # 设置颜色材质
            glEnable(GL_COLOR_MATERIAL)
            glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE)

            # 优化1: 增强光源系统
            # 主光源（提供主要照明和高光）
            light0_position = [5.0, 5.0, 10.0, 1.0]  # 从右前方照射
            light0_ambient = [0.2, 0.2, 0.2, 1.0]
            light0_diffuse = [0.9, 0.9, 0.9, 1.0]
            light0_specular = [1.0, 1.0, 1.0, 1.0]  # 强高光

            glLightfv(GL_LIGHT0, GL_POSITION, light0_position)
            glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient)
            glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse)
            glLightfv(GL_LIGHT0, GL_SPECULAR, light0_specular)

            # 填充光源（减少阴影，从左后方照射）
            light1_position = [-5.0, -3.0, 5.0, 1.0]
            light1_ambient = [0.1, 0.1, 0.1, 1.0]
            light1_diffuse = [0.4, 0.4, 0.4, 1.0]
            light1_specular = [0.2, 0.2, 0.2, 1.0]

            glLightfv(GL_LIGHT1, GL_POSITION, light1_position)
            glLightfv(GL_LIGHT1, GL_AMBIENT, light1_ambient)
            glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse)
            glLightfv(GL_LIGHT1, GL_SPECULAR, light1_specular)

            # 优化2: 添加雾效增强深度感
            glEnable(GL_FOG)
            fog_color = [0.85, 0.85, 0.85, 1.0]  # 与背景色匹配
            glFogfv(GL_FOG_COLOR, fog_color)
            glFogi(GL_FOG_MODE, GL_LINEAR)
            glFogf(GL_FOG_START, 1.0)
            glFogf(GL_FOG_END, 15.0)
            glFogf(GL_FOG_DENSITY, 0.35)

            # 启用反走样
            glEnable(GL_LINE_SMOOTH)
            glEnable(GL_BLEND)
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
            glHint(GL_LINE_SMOOTH_HINT, GL_NICEST)

            # 优化3: 设置材质高光属性
            glMaterialfv(GL_FRONT, GL_SPECULAR, [1.0, 1.0, 1.0, 1.0])
            glMaterialf(GL_FRONT, GL_SHININESS, 80.0)  # 高光强度

            # 新增：创建圆柱体显示列表
            self.create_bond_cylinders()

            self.initialized = True
            logger.info("OpenGL initialized successfully")

        except Exception as e:
            logger.exception(f"Error initializing OpenGL: {str(e)}")
            QMessageBox.critical(self, "OpenGL Error", f"无法初始化OpenGL: {str(e)}")

    def create_bond_cylinders(self):
        """创建用于绘制化学键的圆柱体显示列表"""
        # 单键
        self.bond_cylinders['single'] = glGenLists(1)
        glNewList(self.bond_cylinders['single'], GL_COMPILE)
        self.draw_cylinder(0.05, 1.0, 20)
        glEndList()

        # 双键
        self.bond_cylinders['double'] = glGenLists(1)
        glNewList(self.bond_cylinders['double'], GL_COMPILE)
        # 两个平行圆柱体
        glPushMatrix()
        glTranslatef(0.1, 0, 0)
        self.draw_cylinder(0.05, 1.0, 20)
        glPopMatrix()

        glPushMatrix()
        glTranslatef(-0.1, 0, 0)
        self.draw_cylinder(0.05, 1.0, 20)
        glPopMatrix()
        glEndList()

        # 三键
        self.bond_cylinders['triple'] = glGenLists(1)
        glNewList(self.bond_cylinders['triple'], GL_COMPILE)
        # 三个平行圆柱体
        glPushMatrix()
        glTranslatef(0.15, 0, 0)
        self.draw_cylinder(0.05, 1.0, 20)
        glPopMatrix()

        glPushMatrix()
        glTranslatef(-0.15, 0, 0)
        self.draw_cylinder(0.05, 1.0, 20)
        glPopMatrix()

        self.draw_cylinder(0.05, 1.0, 20)
        glEndList()

    def draw_cylinder(self, radius, height, slices):
        """绘制圆柱体"""
        quad = gluNewQuadric()
        gluQuadricDrawStyle(quad, GLU_FILL)
        gluQuadricNormals(quad, GLU_SMOOTH)

        # 侧面
        gluCylinder(quad, radius, radius, height, slices, 1)

        # 顶部圆盘
        glPushMatrix()
        glTranslatef(0, 0, height)
        gluDisk(quad, 0, radius, slices, 1)
        glPopMatrix()

        # 底部圆盘
        gluDisk(quad, 0, radius, slices, 1)

    def resizeGL(self, width, height):
        try:
            if not self.initialized:
                return

            glViewport(0, 0, width, height)
            glMatrixMode(GL_PROJECTION)
            glLoadIdentity()
            gluPerspective(45.0, width / height, 0.1, 100.0)
            glMatrixMode(GL_MODELVIEW)
        except Exception as e:
            logger.exception(f"Error in resizeGL: {str(e)}")

    def paintGL(self):
        try:
            if not self.initialized:
                return

            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
            glLoadIdentity()

            # 应用平移和旋转
            glTranslatef(self.xTrans, self.yTrans, self.zTrans)
            glRotatef(self.xRot / 16.0, 1.0, 0.0, 0.0)
            glRotatef(self.yRot / 16.0, 0.0, 1.0, 0.0)
            glRotatef(self.zRot / 16.0, 0.0, 0.0, 1.0)

            # 绘制分子
            self.draw_molecule()

            # 在绘制完3D场景后，计算原子屏幕位置（用于标签渲染）
            self.calculate_atom_screen_positions()

            # 在3D场景上叠加2D文本标签
            if self.show_atom_labels and self.atoms:
                # 修复1: 保存OpenGL状态
                glPushAttrib(GL_ALL_ATTRIB_BITS)
                glDisable(GL_DEPTH_TEST)
                glDisable(GL_LIGHTING)
                glMatrixMode(GL_PROJECTION)
                glPushMatrix()
                glLoadIdentity()
                gluOrtho2D(0, self.width(), 0, self.height())
                glMatrixMode(GL_MODELVIEW)
                glPushMatrix()
                glLoadIdentity()

                painter = QPainter(self)
                painter.setRenderHint(QPainter.TextAntialiasing)
                painter.setRenderHint(QPainter.SmoothPixmapTransform)

                # 设置字体
                font = QFont("Arial", 10, QFont.Bold)
                painter.setFont(font)

                for i, (screen_x, screen_y) in enumerate(self.atom_screen_positions):
                    atom_type = self.atom_types[i]

                    # 跳过氢原子（如果设置）
                    if atom_type == 'H' and not self.show_hydrogen_labels:
                        continue

                    # 修复2: 改进文本显示
                    # 设置文本颜色（与原子颜色对比）
                    color = ATOM_COLORS.get(atom_type, (0.5, 0.5, 0.5))
                    # 计算亮度，决定使用黑色还是白色文本
                    brightness = 0.299 * color[0] + 0.587 * color[1] + 0.114 * color[2]
                    text_color = Qt.black if brightness > 0.5 else Qt.white

                    # 计算文本尺寸
                    text_rect = painter.fontMetrics().boundingRect(atom_type)

                    # 添加背景矩形（半透明）
                    bg_color = QColor(0, 0, 0, 180)  # 半透明黑色
                    if brightness > 0.5:  # 如果原子颜色较亮
                        bg_color = QColor(255, 255, 255, 180)  # 半透明白色

                    # 增加偏移量，避免与原子球体重叠
                    label_x = int(screen_x) + 10
                    label_y = int(screen_y) - 10

                    # 绘制背景矩形
                    painter.setPen(Qt.NoPen)
                    painter.setBrush(QBrush(bg_color))
                    painter.drawRect(label_x - 2, label_y - text_rect.height() + 2,
                                     text_rect.width() + 4, text_rect.height() + 2)

                    # 绘制文本
                    painter.setPen(text_color)
                    painter.drawText(label_x, label_y, atom_type)

                painter.end()

                # 修复3: 恢复OpenGL状态
                glMatrixMode(GL_PROJECTION)
                glPopMatrix()
                glMatrixMode(GL_MODELVIEW)
                glPopMatrix()
                glPopAttrib()
        except Exception as e:
            logger.exception(f"Error in paintGL: {str(e)}")

    def calculate_atom_screen_positions(self):
        """计算所有原子在屏幕上的2D位置"""
        self.atom_screen_positions = []

        if not self.atoms:
            return

        # 获取当前OpenGL矩阵
        modelview = glGetDoublev(GL_MODELVIEW_MATRIX)
        projection = glGetDoublev(GL_PROJECTION_MATRIX)
        viewport = glGetIntegerv(GL_VIEWPORT)

        # 为每个原子计算屏幕位置
        for atom in self.atoms:
            x, y, z = atom
            # 使用gluProject将3D坐标转换为2D屏幕坐标
            win_x, win_y, win_z = gluProject(x, y, z, modelview, projection, viewport)
            # OpenGL的Y坐标是从底部开始的，需要转换为Qt的Y坐标（从顶部开始）
            screen_y = self.height() - win_y
            self.atom_screen_positions.append((win_x, screen_y))

    def draw_molecule(self):
        try:
            # 绘制键
            for bond in self.bonds:
                i, j = bond
                x1, y1, z1 = self.atoms[i]
                x2, y2, z2 = self.atoms[j]

                # 计算方向向量
                dx, dy, dz = x2 - x1, y2 - y1, z2 - z1
                length = np.sqrt(dx * dx + dy * dy + dz * dz)

                # 计算旋转角度
                angle = np.degrees(np.arctan2(np.sqrt(dx * dx + dy * dy), dz))

                # 修复1: 明确指定axis为浮点类型数组
                axis = np.array([-dy, dx, 0], dtype=np.float64)
                axis_length = np.sqrt(np.sum(axis * axis))

                # 修复2: 特殊处理Z轴方向的键
                if axis_length < 1e-6:  # 几乎为零，表示键沿Z轴
                    if dz > 0:
                        axis = np.array([1.0, 0.0, 0.0], dtype=np.float64)  # 使用X轴作为旋转轴
                    else:
                        axis = np.array([0.0, 1.0, 0.0], dtype=np.float64)  # 使用Y轴作为旋转轴
                    axis_length = 1.0

                if axis_length > 0:
                    axis /= axis_length  # 现在可以安全地进行浮点除法

                # 绘制键
                glPushMatrix()
                glTranslatef(x1, y1, z1)

                if axis_length > 0:
                    glRotatef(angle, axis[0], axis[1], axis[2])

                # 缩放圆柱体长度
                glScalef(1, 1, length)

                # 选择键类型
                bond_key = 'single'  # 简化处理，所有键都用单键
                if bond_key in self.bond_cylinders:
                    glCallList(self.bond_cylinders[bond_key])
                else:
                    glCallList(self.bond_cylinders['single'])

                glPopMatrix()

            # 绘制原子球体
            for i, atom in enumerate(self.atoms):
                x, y, z = atom
                atom_type = self.atom_types[i]

                # 获取原子颜色
                color = ATOM_COLORS.get(atom_type, (0.5, 0.5, 0.5))
                glColor3f(*color)

                # 创建球体
                quad = gluNewQuadric()
                gluQuadricDrawStyle(quad, GLU_FILL)
                gluQuadricNormals(quad, GLU_SMOOTH)

                glPushMatrix()
                glTranslatef(x, y, z)

                # 根据原子类型调整半径
                if atom_type == 'H':
                    radius = 0.25
                    slices = 16
                    stacks = 16
                else:
                    radius = 0.35
                    slices = 24
                    stacks = 24

                # 优化：增加深度提示
                modelview = glGetDoublev(GL_MODELVIEW_MATRIX)
                z_view = modelview[2][0] * x + modelview[2][1] * y + modelview[2][2] * z + modelview[2][3]
                size_factor = max(0.8, 1.0 - 0.03 * abs(z_view))
                radius *= size_factor

                gluSphere(quad, radius, slices, stacks)
                glPopMatrix()
        except Exception as e:
            logger.exception(f"Error drawing molecule: {str(e)}")

    def mousePressEvent(self, event):
        self.lastX = event.x()
        self.lastY = event.y()

    def mouseMoveEvent(self, event):
        try:
            dx = event.x() - self.lastX
            dy = event.y() - self.lastY

            if event.buttons() & Qt.LeftButton:
                # 旋转
                self.xRot += 8 * dy
                self.yRot += 8 * dx
                self.update()
            elif event.buttons() & Qt.RightButton:
                # 平移
                self.xTrans += dx / 100.0
                self.yTrans -= dy / 100.0
                self.update()

            self.lastX = event.x()
            self.lastY = event.y()
        except Exception as e:
            logger.exception(f"Error in mouseMoveEvent: {str(e)}")

    def wheelEvent(self, event):
        try:
            # 缩放
            num_degrees = event.angleDelta().y() / 8
            num_steps = num_degrees / 15
            self.zTrans += num_steps
            self.update()
        except Exception as e:
            logger.exception(f"Error in wheelEvent: {str(e)}")

    # 新增方法：设置分子数据
    def set_molecule_data(self, atoms, bonds, atom_types):
        """设置要显示的分子数据"""
        self.atoms = atoms
        self.bonds = bonds
        self.atom_types = atom_types

        # 重置视图
        self.xRot = 0
        self.yRot = 0
        self.zRot = 0
        self.xTrans = 0
        self.yTrans = 0
        self.zTrans = -5

        self.update()


class MoleculeApp(QMainWindow):
    def open_in_molview(self):
        """在浏览器中使用 MolView 打开当前分子"""
        try:
            # 获取当前输入内容
            user_input = self.molecule_input.text().strip()
            if not user_input:
                QMessageBox.warning(self, "输入为空", "请输入分子名称或SMILES")
                return

            # 判断是直接使用输入内容作为 SMILES 还是查数据库
            is_smiles_mode = self.smiles_mode.isChecked()
            smiles_to_use = user_input if is_smiles_mode else MOLECULE_DB.get(user_input.lower(), None)

            if not smiles_to_use:
                # 尝试从名称模糊匹配
                for name, smi in MOLECULE_DB.items():
                    if user_input.lower() in name or name in user_input.lower():
                        smiles_to_use = smi
                        break

            if not smiles_to_use:
                reply = QMessageBox.question(
                    self, "未知分子",
                    f"无法识别 '{user_input}'，是否直接以该字符串作为SMILES发送到MolView？",
                    QMessageBox.Yes | QMessageBox.No
                )
                if reply == QMessageBox.Yes:
                    smiles_to_use = user_input
                else:
                    return

            # 编码并生成 URL
            import urllib.parse
            encoded_smiles = urllib.parse.quote(smiles_to_use)
            url = f"https://molview.org/?smiles={encoded_smiles}"

            # 打开浏览器
            webbrowser.open(url)

            # 提示用户
            QMessageBox.information(
                self, "已打开",
                f"正在浏览器中打开 MolView...\n"
                f"SMILES: {smiles_to_use}\n"
                f"网址: {url}"
            )

        except Exception as e:
            QMessageBox.critical(self, "错误", f"无法打开MolView: {str(e)}")

    def __init__(self):
        super().__init__()

        self.setWindowTitle("化学分子结构可视化")
        self.setGeometry(100, 100, 1000, 750)

        # 创建主部件
        main_widget = QWidget()
        self.setCentralWidget(main_widget)

        # 创建主布局
        main_layout = QHBoxLayout(main_widget)

        # 左侧控制面板
        control_panel = QWidget()
        control_panel.setFixedWidth(300)
        control_layout = QVBoxLayout(control_panel)

        # 标题
        title_label = QLabel("化学分子结构可视化")
        title_label.setFont(QFont("Arial", 14))  # 14号字体刚刚好
        title_label.setStyleSheet("color: #2c3e50;")
        title_label.setAlignment(Qt.AlignCenter)
        control_layout.addWidget(title_label)

        # 添加间隔
        control_layout.addSpacing(20)

        # 输入模式选择
        input_mode_label = QLabel("输入模式:")
        control_layout.addWidget(input_mode_label)

        # 创建按钮组
        self.input_mode_group = QButtonGroup(self)

        # 名称输入模式
        self.name_mode = QRadioButton("分子名称")
        self.name_mode.setChecked(True)
        self.input_mode_group.addButton(self.name_mode, 0)

        # SMILES输入模式
        self.smiles_mode = QRadioButton("SMILES字符串")
        self.input_mode_group.addButton(self.smiles_mode, 1)

        mode_layout = QHBoxLayout()
        mode_layout.addWidget(self.name_mode)
        mode_layout.addWidget(self.smiles_mode)
        control_layout.addLayout(mode_layout)

        # 分子名称输入
        input_layout = QHBoxLayout()
        self.molecule_input = QLineEdit()
        self.molecule_input.setPlaceholderText("输入分子名称 (如: 水, 甲烷)")
        input_layout.addWidget(self.molecule_input)

        # 创建按钮
        self.generate_btn = QPushButton("生成")
        self.generate_btn.setStyleSheet("""
            QPushButton {
                background-color: #3498db;
                color: white;
                border: none;
                padding: 5px;
                border-radius: 3px;
            }
            QPushButton:hover {
                background-color: #2980b9;
            }
        """)
        self.generate_btn.clicked.connect(self.generate_molecule)
        input_layout.addWidget(self.generate_btn)

        control_layout.addLayout(input_layout)

        # 添加常用分子选择
        control_layout.addWidget(QLabel("常用分子:"))
        self.molecule_combo = QComboBox()
        self.molecule_combo.addItems([
            "水", "甲烷", "乙烷", "乙烯", "乙炔",
            "二氧化碳", "苯", "氨", "甲醇", "乙醇", "葡萄糖",
            "六氟化硫"  # 新增：添加六氟化硫到常用分子列表
        ])
        control_layout.addWidget(self.molecule_combo)

        # 添加"使用所选"按钮
        self.use_selected_btn = QPushButton("使用所选分子")
        self.use_selected_btn.setStyleSheet("""
            QPushButton {
                background-color: #27ae60;
                color: white;
                border: none;
                padding: 5px;
                border-radius: 3px;
            }
            QPushButton:hover {
                background-color: #219653;
            }
        """)
        self.use_selected_btn.clicked.connect(self.use_selected_molecule)
        control_layout.addWidget(self.use_selected_btn)
        # ==== 新增：在MolView中查看 ====
        self.molview_btn = QPushButton("🌐 在MolView中查看")
        self.molview_btn.setStyleSheet("""
            QPushButton {
                background-color: #16a085;
                color: white;
                border: none;
                padding: 6px;
                border-radius: 3px;
            }
            QPushButton:hover {
                background-color: #13846d;
            }
        """)
        self.molview_btn.clicked.connect(self.open_in_molview)
        control_layout.addWidget(self.molview_btn)

        # 新增：显示选项
        control_layout.addWidget(QLabel("显示选项:"))

        self.show_labels_checkbox = QCheckBox("显示原子标签")
        self.show_labels_checkbox.setChecked(True)
        self.show_labels_checkbox.stateChanged.connect(self.toggle_atom_labels)
        control_layout.addWidget(self.show_labels_checkbox)

        self.show_hydrogen_labels_checkbox = QCheckBox("显示氢原子标签")
        self.show_hydrogen_labels_checkbox.setChecked(False)
        self.show_hydrogen_labels_checkbox.stateChanged.connect(self.toggle_hydrogen_labels)
        control_layout.addWidget(self.show_hydrogen_labels_checkbox)
        # 新增结束

        # 添加间隔
        control_layout.addSpacing(20)

        # 信息面板
        info_frame = QFrame()
        info_frame.setFrameShape(QFrame.StyledPanel)
        info_frame.setStyleSheet("background-color: #f8f9fa; border-radius: 5px;")
        info_layout = QVBoxLayout(info_frame)

        info_label = QLabel("使用说明:")
        info_label.setFont(QFont("Arial", 10, QFont.Bold))
        info_label.setStyleSheet("color: #2c3e50;")
        info_layout.addWidget(info_label)

        instructions = QLabel(
            "• 左键拖动: 旋转分子\n"
            "• 右键拖动: 平移分子\n"
            "• 滚轮: 缩放\n"
            "• 输入分子名称或SMILES字符串"
        )
        instructions.setWordWrap(True)
        instructions.setStyleSheet("background-color: transparent;")
        info_layout.addWidget(instructions)

        control_layout.addWidget(info_frame)

        # 添加分子信息
        control_layout.addWidget(QLabel("当前分子信息:"))
        self.info_text = QLabel("正在初始化...")
        self.info_text.setWordWrap(True)
        self.info_text.setStyleSheet("background-color: #f8f9fa; padding: 5px; border-radius: 3px;")
        control_layout.addWidget(self.info_text)

        # 分子属性标签
        control_layout.addWidget(QLabel("分子属性:"))

        # 创建属性表格
        self.properties_table = QTableWidget(10, 2)
        self.properties_table.setHorizontalHeaderLabels(["属性", "值"])
        self.properties_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.properties_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        self.properties_table.setEditTriggers(QTableWidget.NoEditTriggers)  # 禁止编辑
        self.properties_table.setSelectionBehavior(QTableWidget.SelectRows)  # 选择整行
        self.properties_table.setStyleSheet("""
            QTableWidget {
                background-color: #f8f9fa;
                border-radius: 3px;
                gridline-color: #e0e0e0;
            }
            QTableWidget::item {
                padding: 4px;
            }
        """)
        control_layout.addWidget(self.properties_table, 1)

        # 添加伸缩空间
        control_layout.addStretch()

        # 右侧分子视图
        self.viewer = MoleculeViewer()

        # 将部件添加到主布局
        main_layout.addWidget(control_panel)
        main_layout.addWidget(self.viewer, 1)

        # 使用定时器确保OpenGL上下文初始化完成
        self.init_timer = QTimer()
        self.init_timer.timeout.connect(self.initialize_viewer)
        self.init_timer.start(100)

    def initialize_viewer(self):
        """初始化分子视图"""
        if self.viewer.initialized:
            self.init_timer.stop()
            # 生成默认分子（水）
            self.generate_molecule("水")
        else:
            logger.info("Waiting for OpenGL initialization...")

    # 新增方法：创建完美的正八面体结构
    def create_octahedral_structure(self, center_atom, ligand_atom, bond_length=1.56):
        """
        创建完美的正八面体分子结构

        参数:
        center_atom - 中心原子符号 (如 "S")
        ligand_atom - 配体原子符号 (如 "F")
        bond_length - 键长 (Å)，默认为S-F键长1.56Å

        返回:
        (success, mol, atoms, bonds, atom_types)
        """
        try:
            logger.info(f"Creating perfect octahedral structure for {center_atom}{ligand_atom}6")

            # 正八面体顶点坐标 (中心在原点，顶点在坐标轴上)
            # 使用精确计算确保所有键等长
            positions = [
                (bond_length, 0, 0),  # +x
                (-bond_length, 0, 0),  # -x
                (0, bond_length, 0),  # +y
                (0, -bond_length, 0),  # -y
                (0, 0, bond_length),  # +z
                (0, 0, -bond_length)  # -z
            ]

            # 创建原子列表
            atoms = [(0.0, 0.0, 0.0)]  # 中心原子
            atoms.extend(positions)  # 6个配位原子

            # 创建原子类型列表
            atom_types = [center_atom] + [ligand_atom] * 6

            # 创建化学键 (所有键都是中心原子到配位原子)
            bonds = [(0, i) for i in range(1, 7)]

            # 创建RDKit分子对象用于属性计算
            from rdkit import Chem
            mol = Chem.RWMol()

            # 添加中心原子
            center_idx = mol.AddAtom(Chem.Atom(center_atom))

            # 添加配位原子和键
            for i in range(6):
                ligand_idx = mol.AddAtom(Chem.Atom(ligand_atom))
                mol.AddBond(center_idx, ligand_idx, Chem.BondType.SINGLE)

            mol = mol.GetMol()

            # ===== 关键修复：确保分子规范化 =====
            try:
                Chem.SanitizeMol(mol)
            except:
                pass  # 如果规范化失败，继续

            # 创建3D构型
            conf = Chem.Conformer(mol.GetNumAtoms())
            conf.SetAtomPosition(0, (0.0, 0.0, 0.0))  # 中心原子

            for i, pos in enumerate(positions):
                conf.SetAtomPosition(i + 1, pos)

            mol.AddConformer(conf)

            return True, mol, atoms, bonds, atom_types

        except Exception as e:
            logger.exception(f"Error creating octahedral structure: {str(e)}")
            return False, None, [], [], []

    # 修改 generate_molecule 方法
    def generate_molecule(self, name=None):
        # 确保name是字符串或None
        if name is not None and not isinstance(name, str):
            logger.error(f"generate_molecule called with non-string name: {type(name)}")
            name = None

        if name is None:
            name = self.molecule_input.text().strip()
            if not name:
                QMessageBox.warning(self, "输入错误", "请输入分子名称或SMILES字符串")
                return

        # 检查输入模式
        is_smiles = self.smiles_mode.isChecked()

        # 关键修改：直接在MoleculeApp中生成分子数据
        success, mol, atoms, bonds, atom_types = self.generate_molecule_data(name, is_smiles)
        if success:
            # 将分子数据传递给viewer
            self.viewer.set_molecule_data(atoms, bonds, atom_types)

            # 更新基本信息
            input_type = "SMILES" if is_smiles else "名称"
            self.info_text.setText(f"已加载: {name} ({input_type})\n原子数: {len(atoms)}")
            logger.info(f"Successfully loaded molecule: {name}")

            # 计算并显示分子属性
            self.display_molecular_properties(mol, name)
        else:
            error_msg = f"无法生成分子: {name}"
            if not is_smiles:
                error_msg += "\n请尝试:\n- 检查拼写\n- 使用标准名称\n- 切换到SMILES模式直接输入结构式"
            else:
                error_msg += "\n请检查SMILES格式是否正确"

            self.info_text.setText(f"加载失败: {name}\n请检查输入")
            logger.warning(f"Failed to load molecule: {name}")
            QMessageBox.warning(self, "分子生成失败", error_msg)

    # 新增方法：在MoleculeApp中生成分子数据
    def generate_molecule_data(self, name, is_smiles=False):
        """在MoleculeApp中生成分子数据，而不是在MoleculeViewer中"""
        try:
            # 确保name是字符串
            if not isinstance(name, str):
                logger.error(f"Expected name to be a string, but got {type(name)}: {name}")
                return False, None, [], [], []

            logger.info(f"Generating molecule: {name}, is_smiles={is_smiles}")

            # 尝试导入RDKit
            try:
                from rdkit import Chem
                from rdkit.Chem import AllChem
            except Exception as e:
                logger.exception(f"Failed to import RDKit: {str(e)}")
                QMessageBox.critical(self, "RDKit错误",
                                     "无法加载RDKit库。\n"
                                     "这通常是NumPy版本不兼容导致的。\n"
                                     "请确保安装了NumPy 1.x版本:\n"
                                     "运行: pip install 'numpy<2'")
                return False, None, [], [], []

            # 检查是否是请求正八面体结构
            name_lower = name.lower().strip()
            if "八面体" in name_lower or "octahedral" in name_lower or name_lower in ["sf6", "六氟化硫",
                                                                                      "sulfur hexafluoride"]:
                # 专门创建完美的正八面体结构
                if "硫" in name_lower or "sulfur" in name_lower or name_lower in ["sf6", "六氟化硫",
                                                                                  "sulfur hexafluoride"]:
                    return self.create_octahedral_structure("S", "F", 1.56)
                elif "钼" in name_lower or "molybdenum" in name_lower:
                    return self.create_octahedral_structure("Mo", "CO", 2.1)  # Mo(CO)6
                elif "铬" in name_lower or "chromium" in name_lower:
                    return self.create_octahedral_structure("Cr", "Cl", 2.3)  # CrCl6^3-
                else:
                    # 默认创建硫-氟正八面体
                    return self.create_octahedral_structure("S", "F", 1.56)

            # 如果是SMILES直接使用
            if is_smiles:
                smiles = name
                logger.info(f"Using SMILES directly: {smiles}")
            else:
                # 尝试匹配数据库
                if name_lower in MOLECULE_DB:
                    smiles = MOLECULE_DB[name_lower]
                    logger.info(f"Found in database: {name} -> {smiles}")
                else:
                    # 尝试作为SMILES处理
                    smiles = name
                    logger.info(f"Not found in database, trying as SMILES: {smiles}")

            # 从SMILES创建分子
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                # ===== 关键修复：确保分子规范化 =====
                try:
                    Chem.SanitizeMol(mol)
                except:
                    pass  # 如果规范化失败，继续

                # 添加氢原子
                mol = Chem.AddHs(mol)

                # 生成3D构型
                try:
                    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
                    AllChem.UFFOptimizeMolecule(mol)
                except:
                    # 如果3D生成失败，尝试2D
                    AllChem.Compute2DCoords(mol)

                # 提取原子和坐标
                atoms = []
                atom_types = []
                for atom in mol.GetAtoms():
                    pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
                    atoms.append((pos.x, pos.y, pos.z))
                    atom_types.append(atom.GetSymbol())

                # 提取化学键
                bonds = []
                for bond in mol.GetBonds():
                    bonds.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))

                logger.info(f"Molecule generated successfully: {name}")
                return True, mol, atoms, bonds, atom_types
            else:
                logger.error(f"Failed to create molecule from SMILES: {smiles}")
                # 尝试作为名称再次匹配（更宽松的匹配）
                for key in MOLECULE_DB:
                    if name_lower in key.lower() or key.lower() in name_lower:
                        return self.generate_molecule_data(key, is_smiles=False)

                return False, None, [], [], []
        except Exception as e:
            logger.exception(f"Error generating molecule {name}: {str(e)}")
            return False, None, [], [], []

    def use_selected_molecule(self):
        name = self.molecule_combo.currentText()
        self.molecule_input.setText(name)
        self.generate_molecule()

    # 新增方法：控制标签显示
    def toggle_atom_labels(self, state):
        """切换是否显示原子标签"""
        self.viewer.show_atom_labels = bool(state)
        self.viewer.update()

    def toggle_hydrogen_labels(self, state):
        """切换是否显示氢原子标签"""
        self.viewer.show_hydrogen_labels = bool(state)
        self.viewer.update()

    def display_molecular_properties(self, mol, name):
        """显示分子属性"""
        try:
            # 尝试导入RDKit
            from rdkit import Chem
            from rdkit.Chem import AllChem, Descriptors, Lipinski, rdMolDescriptors

            properties = calculate_molecular_properties(mol)

            # 设置表格行数
            self.properties_table.setRowCount(len(properties))

            # 填充表格
            row = 0
            for prop, value in properties.items():
                prop_item = QTableWidgetItem(prop)
                value_item = QTableWidgetItem(str(value))

                # 设置样式
                prop_item.setBackground(QBrush(QColor(240, 240, 240)))
                prop_item.setFont(QFont("Arial", 9, QFont.Bold))

                self.properties_table.setItem(row, 0, prop_item)
                self.properties_table.setItem(row, 1, value_item)
                row += 1

            # 更新标题
            input_type = "SMILES" if self.smiles_mode.isChecked() else "名称"
            self.info_text.setText(f"已加载: {name} ({input_type})\n原子数: {len(self.viewer.atoms)}")
        except Exception as e:
            logger.exception(f"Error displaying molecular properties: {str(e)}")
            self.info_text.setText(f"加载分子: {name}\n属性计算失败\n请确保安装了兼容的NumPy版本")


if __name__ == "__main__":
    # 检查依赖
    missing_deps = check_dependencies()
    if missing_deps:
        error_msg = "缺少必要的依赖库:\n\n" + "\n".join(missing_deps) + "\n\n请先安装这些库。"
        logger.error(error_msg)
        app = QApplication(sys.argv)
        QMessageBox.critical(None, "依赖错误", error_msg)
        sys.exit(1)

    try:
        app = QApplication(sys.argv)
        window = MoleculeApp()
        window.show()
        logger.info("Application started successfully")
        sys.exit(app.exec_())
    except Exception as e:
        logger.exception("Application crashed")
        QMessageBox.critical(None, "应用程序错误", f"程序发生严重错误: {str(e)}")
        sys.exit(1)
