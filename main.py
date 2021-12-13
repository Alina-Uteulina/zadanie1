import numpy as np

p_L = int(input())
lambda_L = int(input())
p_g = int(input())
m_L = int(input())
m_g = int(input())
g = 9.8
teta = 90
f_tr = int(input())
v_s = int(input())
d = int(input())
q_L = int(input())
A_p = int(input())
#параметры двухфазного потока
p_tr = p_L*lambda_L + p_g*(1-lambda_L)
m_tr = m_L*lambda_L + m_g*(1-lambda_L)
v_sL = q_L/A_p
v_sg = 0.25*v_s + 0.333*v_sL
v_tr = v_sL + v_sg
#градиент давления для пузырькового режима 
class puz(p_tr,g,f_tr,v_tr,d,m_tr):    
    def _init_(puz):
        grad_puz = (p_tr*g*np.sin(teta))+(f_tr*p_tr*v_tr**2/2*d)
    return g
puz(p_tr,g,f_tr,v_tr,d,m_tr)
#Пробковый режим
beta = int(input())
H_LLs = int(input())
f_Ls = int(input())
m_Ls = int(input())
v_gtb = int(input())
v_gLs = int(input())
v_Ltb = int(input())
H_LLs = int(input())
#параметры
p_Ls = p_L*H_LLs + p_g*(1-H_LLs)
v_sg1 = beta*v_gtb*(1-H_Ltb)+(1-beta)*v_gLs*(1-H_LLs)
v_Ls1 = (1-beta)*v_LLs*H_LLs-beta*v_Ltb*H_Ltb
v_m = v_sg1 + v_sL1
#градиент давления для пробкового режима
class prob(f_Ls,p_Ls1,v_m,d,beta,p_g,teta,g)
    def _init_(prob):
        grad_prob = ((1-beta)*p_Ls1+beta*p_g)*g*np.sin(teta)+f_Ls*p_Ls*v_m**2/2*d*(1-beta)
    return grad_prob
prob(f_Ls,p_Ls1,v_m,d,beta,p_g,teta,g)
#параметры для эмульсионного режима
C_0 = 1.15
p__L = int(input())
lambda__L = int(input())
p__g = int(input())
m__L = int(input())
m__g = int(input())
g = 9.8
teta = 90
f__tr = int(input())
v__s = int(input())
d = int(input())
q__L = int(input())
A__p = int(input())
#параметры двухфазного потока
p__tr = p__L*lambda__L + p__g*(1-lambda__L)
m__tr = m__L*lambda__L + m__g*(1-lambda__L)
v__tr = v__sL + v__sg
v__sL = q__L/A__p
v__sg = np.sin(teta)/(4-C_0)*C_0*v__sL+v__s
#градиент давления для эмульсионного режима
class mus(p__tr,g,f__tr,v__tr,d,teta):
    def _init_(mus):
        grad_mus = (p__tr*g*np.sin(teta))+(f__tr*p__tr*v__tr**2/2*d)
    return grad_mus
mus(p__tr,g,f__tr,v__tr,d,teta)
#Кольцевой режим
p_L3 = int(input())
p_c = int(input())
p_g3 = int(input())
sigma_L = int(input())
m_g = int(input())
g = 9.8
teta = 90
f_sc = int(input())
d = int(input())
q_L = int(input())
q_g = int(input())
A_p = int(input())
delta = 0.1
#параметры двухфазного потока
v_sL3 = q_L/A_p
v_sg3 = q_g/A_p
v_kr = 10000*v_sg3*m_g/sigma_L*(p_g3/p_L3)**0.5
F_e = 1 - mh.exp*(-0.125(v_kr - 1.5))
v_sc = F_e*v_sL3 + v_sg3
dp = f_sc*p_c*v_sc**2/2*d
if F_e > 0.9:
    z = 1+300*delta
else:
    z = 1 + 24*delta*(p_L/p_g)**(1/3)
dp_c = z/(1-2*delta)**5*dp+p_c*g*np.sin(teta)
fi = (dp_c)-p_c*g*np.sin(teta)/dp
#градиент давления для кольцевого режима
class kol(fi,dp,g,p_c,teta):
    def _iniy_(kol):
        grad_kol = fi*dp+g*p_c*np.sin(teta)
    return grad_kol
kol(fi,dp,g,p_c,teta)