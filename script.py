

import tellurium as te
import sys
import os
from pathlib import Path
dir_file = Path(__file__).resolve().parent
main_dir = os.path.join(dir_file)
sys.path.insert(0,main_dir)
import json
sys.path.insert(0,main_dir)
from data.observations import observations,t2m,select_obs,packages
from models.params import fixed_params
from tools import dirs, activation
# from models.models import Macrophage
import matplotlib
import matplotlib.pyplot as plt
import tellurium as te
from models.models import Macrophage
params ={
    "IL8": 4.724288120656182,
    "IL8_m": 1.2138313631672166,
    "k_il8_p": 0.024397471013379857,
    "k_il8m_il8": 0.14596703636802016,
    "IL8_R": 2.310405727945194,
    "pIL8_R": 6.666584229891921,
    "k_il8r_b": 0.6167265055192227,
    "k_il8r_ub": 0.19347421128048325,
    "k_il8r_a": 0.014547178907135416,
    "k_il8r_da": 0.46775013717966707,
    "kd_il8_irak_p": 296092.5670587952,
    "k_il8_irak_p": 751303.2754822965,
    "o_il8_irak_p": 0.205898666571927,
    "k_nfkb_il8_p": 1.2380664811373663,
    "kd_nfkb_il8_p": 768.6185511351218,
    "o_nfkb_il8_p": 0.8408046136068019,
    "k_ap1_il8_p": 1.2774381921529425,
    "kd_ap1_il8_p": 786.3889374437415,
    "o_ap1_il8_p": 0.9967236511541081,
    "k_rho_a": 1.1167951937654834,
    "kd_rho_a": 2456.098246911869,
    "k_rho_da": 0.6171084280662366,
    "k_rho_nfkb_a": 22697.33985748905,
    "kd_rho_nfkb_a": 42339.85375059342,
    "o_rho_nfkb_a": 0.0024194362951256987,
    "k_rho_pi3k_a": 23821.182456496528,
    "kd_rho_pi3k_a": 22183.11093689271,
    "o_rho_pi3k_a": 0.2409415621442279,
    "k_rho_stat3_a": 5743.654088652801,
    "kd_rho_stat3_a": 20174.158326999215,
    "o_rho_stat3_a": 0.9991400308481002,
    "_pSTAT3": 45.232509476466234,
    "k_pstat3_d": 0.5357996204883932,
    "k_pstat3_b": 0.6563465730041598,
    "k_pstat3_ub": 0.5224575788229661,
    "_pPI3K": 95.73487662437961,
    "k_ppi3k_d": 0.06376488201257813,
    "K_ppi3k_a": 0.9101379676030306,
    "K_ppi3k_da": 0.5488009580354664
}
## activate macrophage
if False:
    activation.activation_LPS()
## curate inferred parameters to insert into the model
if False:
    file_name = os.path.join(main_dir,'outputs','inferred_params_ILs_1.json')
    with open(file_name,'r') as f:
        ss = f.read()
    ss = ss.replace(",",";")
    ss = ss.replace(":","=")
    ss = ss.replace("\""," ")
    print(ss)
if True:
    ss = str(params)
    ss = ss.replace(",",";\n")
    ss = ss.replace(":","=")
    ss = ss.replace("\""," ")
    ss = ss.replace("\'"," ")
    print(ss)