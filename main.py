from расчеты import *
from scipy.integrate import solve_ivp

class Gradient:

    def __init__(self, d):
        self.d = None

    @staticmethod
    def calc_fp(v_sl, fp):
        """
        Определение структуры потока
        :param v_sl: скорость жидкости
        :return: номер режима потока, безразмерн.
                режим потока:
                * 1 - пузырьковый;
                * 2 - пробковый;
                * 3 - эмульсионный;
                * 4 - кольцевой;
        """
        sg1 == 0.25 * v_s + 0.333 * v_sL
        sg4 == 3.1 * (g * sigma_L * (p_L - p_g) / p_g ** 2) ** 1 / 4

        if v_sl < sg1:
            fp = 1
        elif v_m > sg1:
            fp = 2
        elif v__sg < sg4:
            fp = 3
        elif v_sc > sg4:
            fp = 4
        return fp


    def puz(p_tr, g, teta, f_tr, v_tr, d, v_sg, v_sl):
        """
        Функция расчета градиента для пузырькового режима

        Parameters
        ----------
        :param p_tr: плотность
        :param g: коэффициент
        :param teta: угол наклона
        :param f_tr: сила трения
        :param v_tr: скорость двухфазного потока
        :param v_sg: скорость газа
        :param v_sl: скорость жидкости
        :param d: коэффициент
        """
        funct_gpuz = (p_tr * g * np.sin(teta)) # гравитационная составляющая
        funct_tpuz = (f_tr * p_tr * v_tr ** 2 / 2 * d) # составляющая по трению
        grad_puz = funct_gpuz + funct_tpuz
        return grad_puz

    def prob(beta, p_Ls, p_g, g, teta, f_Ls, v_m, d, v_sg, C_0,v_sl, v_s):
        """
        расчет градиента давления для пробкового режима

        Parameters
        ----------
        :param beta: соотношение длины
        :param g: коэффициент
        :param teta: угол наклона трубы
        :param f_Ls: сила трения
        :param v_m: скорость смеси
        :param v_sg: скорость газа
        :param v_sl: скорость жидкости
        :param d: коэффициент
        :param C_0: коэффициент
        :param v_s: скорость скольжения
        """
        funct_gpr = ((1 - beta) * p_Ls + beta * p_g) * g * np.sin(teta) # гравитационная составляющая
        funct_tpr =  f_Ls * p_Ls * v_m ** 2 / 2 * d * (1 - beta) # составляющая по трению
        grad_prob = funct_gpr + funct_tpr
        return grad_prob

    def mus(p_tr, g, teta, f_tr, v_tr, d, v_sg, C_1, v_sl, v_s):
        """
        расчет градиенты давления для эмульсионного режима

        Parameters
        ----------
        :param p_tr: плотность
        :param g: коэффициент
        :param teta: угол наклона трубы
        :param f_tr: сила трения
        :param v_tr: скорость двухфазного потока
        :param v_sg: скорость газа
        :param v_sl: скорость жидкости
        :param d: коэффициент
        :param C_1: коэффициент
        :param v_s: скорость скольжения
        """
        funct_gmus = p_tr * g * np.sin(teta) # гравитационная составляющая
        funct_tmus = (f_tr * p_tr * v_tr ** 2 / 2 * d) # составляющая по трению
        grad_mus = funct_gmus + funct_tmus
        return grad_mus

    def kol(fi, dp, g, p_c, teta, v_sg, sigma_L, p_L, p_g,):
        """
        расчет давления для кольцевого режима

        Parameters
        ----------
        :param fi: коэффициент
        :param g: коэффициент
        :param teta: угол наклона трубы
        :param dp: состовляющая градиента давления по трению для газового ядра
        :param p_c: плотность газового ядра
        :param v_sg: скорость газа
        :param sigma_L: поверхностное натяжение
        :param p_L: коэффициент
        :param p_g: коэффициент
        """
        funct_gkol = fi * dp # гравитационная составляющая
        funct_tkol = g * p_c * np.sin(teta) # составляющая по трению
        grad_kol = funct_gkol + funct_tkol
        return grad_kol

def grad(Ansari):
    dp = Ansari.grad()
    return dp

result = solve_ivp(fun_grad)