import math
import pprint
import numpy as np
import scipy as sp

def get_states(hits,window):
    state_dict = {}
    state = "1"*window
    state_dict[state] = 0
    state_idx = 1
    for i in range(1,(2**window)-1):
        state = ""
        bits = i
        for j in range(window):
            state += "1" if bits & 1 else "0"
            bits >>= 1
        if state.count("1") >= hits:
            state_dict[state] = state_idx
            state_idx += 1
    state_dict["out"] =  state_idx
    return state_dict
    
def optimize1d(state_dict, horizon, confidence, getUtilization):
    results = []
    rranges=[(.5, .9999)]
    params, value, _, _ = sp.optimize.brute(test1d, rranges, Ns=1000, args = (state_dict, horizon, confidence, getUtilization), workers=8, full_output=1, finish=None)
    c = test1d(params, state_dict, horizon, confidence, getUtilization, full_output=True)
    if c < confidence:
        value = math.inf
    return value, params
    
def optimize2d(state_dict, horizon, confidence, getUtilization):
    results = []
    rranges=((.5, .9999), (.5, .9999))
    params, value, _, _ = sp.optimize.brute(test2d, rranges, Ns=1000, args = (state_dict, horizon, confidence, getUtilization), workers=8, full_output=1, finish=None)
    c = test2d(params, state_dict, horizon, confidence, getUtilization, full_output=True)
    if c < confidence:
        value = math.inf
    return value, params

def optimize4d(state_dict, horizon, confidence, getUtilization):
    results = []
    rranges=((.5, .9999), (.5, .9999), (.5, .9999), (.5, .9999))
    params, value, _, _ = sp.optimize.brute(test4d, rranges, Ns=25, args = (state_dict, horizon, confidence, getUtilization), workers=8, full_output=1, finish=None)
    return value, params

def test1d (z, state_dict, horizon, confidence, getUtilization, full_output=False):
    percentile = z
    states = len(state_dict)
    transitions = np.zeros((states,states))
    for state,idx in state_dict.items():
        if state == "out":
            transitions[idx][idx] = 1.0
            continue
        next_hit = state[1:] + "1"
        next_miss = state[1:] + "0"
        next_hit_idx = state_dict[next_hit]
        if next_miss in state_dict:
            next_miss_idx = state_dict[next_miss]
        else:
            next_miss_idx = state_dict["out"]
        transitions[idx][next_hit_idx] = percentile
        transitions[idx][next_miss_idx] = 1-percentile
    dtmc = np.matrix(transitions)
    h_step_transitions = dtmc ** horizon
    c = 1.0 - h_step_transitions.item((0,state_dict["out"]))
    if full_output:
        return c
    if c > confidence:
        return getUtilization(percentile)
    else:
        return math.inf

def test2d(z, state_dict, horizon, confidence, getUtilization, full_output=False):
    p1, p2 = z
    if horizon % 2 != 0:
        raise ValueError("H should be a multiple of 2")
    states = len(state_dict)
    transitions1 = np.zeros((states,states))
    transitions2 = np.zeros((states,states))
    for state,idx in state_dict.items():
        if state == "out":
            transitions1[idx][idx] = 1.0
            transitions2[idx][idx] = 1.0
            continue
        next_hit = state[1:] + "1"
        next_miss = state[1:] + "0"
        next_hit_idx = state_dict[next_hit]
        if next_miss in state_dict:
            next_miss_idx = state_dict[next_miss]
        else:
            next_miss_idx = state_dict["out"]
        transitions1[idx][next_hit_idx] = p1
        transitions1[idx][next_miss_idx] = 1-p1
        transitions2[idx][next_hit_idx] = p2
        transitions2[idx][next_miss_idx] = 1-p2
    dtmc2step = np.matrix(transitions1) * np.matrix(transitions2)
    h_step_transitions = dtmc2step ** (int(horizon/2))
    c = 1.0 - h_step_transitions.item((0,state_dict["out"]))
    if full_output:
        return c
    if c > confidence:
        return getUtilization(p1,p2)
    else:
        return getUtilization(.999)

def test4d(z, state_dict, horizon, confidence, getUtilization):
    p1, p2, p3, p4 = z
    if horizon % 4 != 0:
        raise ValueError("H should be a power of 2")
    if p1>= 1.0 or  p2>= 1.0  or  p3>= 1.0  or  p4>= 1.0:
        print(f"{p1},{p2},{p3},{p4}")
        raise ValueError("PROBLEM!")
    states = len(state_dict)
    transitions1 = np.zeros((states,states))
    transitions2 = np.zeros((states,states))
    transitions3 = np.zeros((states,states))
    transitions4 = np.zeros((states,states))
    for state,idx in state_dict.items():
        if state == "out":
            transitions1[idx][idx] = 1.0
            transitions2[idx][idx] = 1.0
            transitions3[idx][idx] = 1.0
            transitions4[idx][idx] = 1.0
            continue
        next_hit = state[1:] + "1"
        next_miss = state[1:] + "0"
        next_hit_idx = state_dict[next_hit]
        if next_miss in state_dict:
            next_miss_idx = state_dict[next_miss]
        else:
            next_miss_idx = state_dict["out"]
        transitions1[idx][next_hit_idx] = p1
        transitions1[idx][next_miss_idx] = 1-p1
        transitions2[idx][next_hit_idx] = p2
        transitions2[idx][next_miss_idx] = 1-p2
        transitions3[idx][next_hit_idx] = p3
        transitions3[idx][next_miss_idx] = 1-p3
        transitions4[idx][next_hit_idx] = p4
        transitions4[idx][next_miss_idx] = 1-p4
    dtmc3step = np.matrix(transitions1) * np.matrix(transitions2) * np.matrix(transitions3) * np.matrix(transitions4)
    h_step_transitions = dtmc3step ** (int(horizon/4))
    c = 1.0 - h_step_transitions.item((0,state_dict["out"]))
    if getUtilization(p1,p2,p3,p4) <= 0:
        raise("leq 0!")
    if c > confidence:
        return getUtilization(p1,p2,p3,p4)
    else:
        return getUtilization(.999)

def getParetoUtilization(*args):
    quantiles = []
    for p in args:
        #Quantile for finite variance pareto, shifted
        if p <1.0:
            quantile = ((1.0-p)**(-1.0/2.1)) - 1.0
            quantiles.append(quantile)
        else:
            quantiles.append(math.inf)
    return np.mean(quantiles)

def getNormalUtilization(*args):
    quantiles = []
    for p in args:
        #Quantile for std normal
        if p <1.0:
            quantile = sp.stats.norm.ppf(p)
            quantiles.append(quantile)
        else:
            quantiles.append(math.inf)
    return np.mean(quantiles)

def getUniformUtilization(*args):
    quantiles = []
    #Unif(0,10) quanitle function (lol)
    for p in args:
        quantiles.append(10.0*p)
    return np.mean(quantiles)






if __name__ == "__main__":
    print("window,hits,utilization,confidence,dof,distribution")
    for confidence in [.99]:
        for window in range(1,7):
            for hits in range(1, window):
                states = get_states(hits,window)
                utilization, percentile =  optimize1d(states, 100, confidence,getParetoUtilization)
                print(f"{window},{hits},{utilization},{confidence},1,pareto")
                utilization, percentile =  optimize1d(states, 100, confidence,getNormalUtilization)
                print(f"{window},{hits},{utilization},{confidence},1,normal")
                utilization, percentile =  optimize1d(states, 100, confidence,getUniformUtilization)
                print(f"{window},{hits},{utilization},{confidence},1,uniform")
    for confidence in [.99]:
        for window in range(1,7):
            for hits in range(1, window):
                states = get_states(hits,window)
                utilization,percentile =  optimize2d(states, 100, confidence,getParetoUtilization)
                print(f"{window},{hits},{utilization},{confidence},2,pareto")
                utilization,percentile =  optimize2d(states, 100, confidence,getNormalUtilization)
                print(f"{window},{hits},{utilization},{confidence},2,normal")
                utilization,percentile =  optimize2d(states, 100, confidence,getUniformUtilization)
                print(f"{window},{hits},{utilization},{confidence},2,uniform")

   # for confidence in [.99]:
   #     for window in range(1,7):
   #         for hits in range(1, window):
   #             states = get_states(hits,window)
   #             utilization,percentile =  optimize4d(states, 100, confidence,getUniformUtilization)
   #             print(f"{window},{hits},{utilization},{confidence},4,uniform")

   # 
