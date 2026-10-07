from cobaya.yaml import yaml_load
info_txt = r"""
likelihood:
#    sn: import_module('likelihood').snlike
    cc: import_module('likelihood').cclike
    bao: import_module('likelihood').baolike
params:      
    H0: 
        prior: {min: 60,  max: 80}
        latex: H_0
    q0: 
        prior: {min: -1, max: 1}
        latex: q_0
    j0: 
        prior: {min: -20, max: 5}
        latex: j_0
    s0: 
        prior: {min: -170, max: 10}
        latex: s_0
force: True
sampler:
    mcmc:
        burn_in: 100
        max_tries: 50000
        max_samples: 50000
        learn_proposal: True

output: ../chains/cb
"""

info = yaml_load(info_txt)
from cobaya.run import run
updated_info, products = run(info)