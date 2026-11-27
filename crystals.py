'''
A simple library with chosen crystals to simplify importing in other scripts
n1: extraordinary index / first extraordinary wave in biaxial crystals
n2: ordinary index / second extraordinary wave in biaxial crystals

All data comes from https://refractiveindex.info/
'''

# 0.20–5.0 µm
sapphire = {"n1": lambda um: (1+1.4313493/(1-(0.0726631/um)**2)+0.65054713/(1-(0.1193242/um)**2)+5.3414021/(1-(18.028251/um)**2))**.5,
            "n2": lambda um: (1+1.5039759/(1-(0.0740288/um)**2)+0.55069141/(1-(0.1216529/um)**2)+6.5927379/(1-(20.072248/um)**2))**.5}

# 0.198–2.05 µm
quartz = {"n1": lambda um: (1+0.28851804+1.09509924/(1-1.02101864e-2/um**2)+1.15662475/(1-100/um**2))**.5,
          "n2": lambda um: (1+0.28604141+1.07044083/(1-1.00585997e-2/um**2)+1.10202242/(1-100/um**2))**.5}

# 0.204-2.172 µm
calcite = {"n1": lambda um: (1+0.35859695+0.82427830/(1-1.06689543e-2/um**2)+0.14429128/(1-120/um**2))**.5,
           "n2": lambda um: (1+0.73358749+0.96464345/(1-1.94325203e-2/um**2)+1.82831454/(1-120/um**2))**.5}

# 0.214–1.53 µm
kdp = {"n1": lambda um: (2.132668+3.2279924*um**2/(um**2-400)+0.008637494/(um**2-0.0122810))**.5,
       "n2": lambda um: (2.259276+13.00522*um**2/(um**2-400)+0.01008956/(um**2-0.0129426))**.5}
