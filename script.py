

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
from tools import dirs, tools
# from models.models import Macrophage
import matplotlib
import matplotlib.pyplot as plt
import tellurium as te
from models.models import Macrophage

## activate macrophage
if False:
    tools.activation_LPS()
## curate inferred parameters to insert into the model
if True:
    file_name = os.path.join(main_dir,'outputs','inferred_params_IL8.json')
    with open(file_name,'r') as f:
        ss = f.read()
    ss = ss.replace(",",";")
    ss = ss.replace(":","=")
    ss = ss.replace("\""," ")
    print(ss)