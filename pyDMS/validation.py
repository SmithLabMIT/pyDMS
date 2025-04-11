
H2S_prop = {"Ch": 150,"kd":10.8,"b":0.15, "Ch_err":6.6, "kd_err":0.6,"b_err":0.0038}
C3H8_prop = {"Ch": 62.8,"kd": 10.5,"b":2.9, "Ch_err":4, "kd_err":0.4,"b_err":0.1}
C3H6_prop = {"Ch": 57.9,"kd":12.7,"b":2.4, "Ch_err":3.3, "kd_err":0.2,"b_err":0.3}
press = np.linspace(10,10,1)
# CASE 1: simple pure DMS : PASS
dms_pure_C3H8 = dms(press, C3H8_prop)
#print(f'p={dms_pure_C3H8[0]},c={dms_pure_C3H8[1]},c_err={dms_pure_C3H8[2]}')

# CASE 2: Virial pure DMS : conditional pass on fugacity
dms_pure_C3H8 = dms(press, C3H8_prop,eos='virial',T=308, gas = 'C3H8')
#print(f'p={dms_pure_C3H8[0]},c={dms_pure_C3H8[1]},c_err={dms_pure_C3H8[2]}')

# CASE 3: PR pure DMS : conditional pass on fugacity
dms_pure_H2S = dms(press, H2S_prop, eos='pr',T=308, gas = 'H2S')
#print(f'p={dms_pure_H2S[0]},c={dms_pure_H2S[1]},c_err={dms_pure_H2S[2]}')

# CASE 4: 2-gas mixed DMS: PASS
dms_mixed_H2S = dms_mixed(press, [H2S_prop,C3H8_prop],x=[0.4,0.6])
#print(f'p={dms_mixed_H2S[0]},c={dms_mixed_H2S[1]},c_err={dms_mixed_H2S[2]}')

# CASE 5: 2-gas mixed Virial & PR DMS : conditional pass on fugacity
dms_mixed_H2S = dms_mixed(press, [H2S_prop,C3H8_prop],x=[0.4,0.6], eos = ['pr','virial'],T=[308,308],gases=['H2S','C3H8'])
#print(f'p={dms_mixed_H2S[0]},c={dms_mixed_H2S[1]},c_err={dms_mixed_H2S[2]}')

# CASE 6: 3-gas mixed DMS: PASS
dms_mixed_H2S = dms_mixed(press, [H2S_prop,C3H8_prop,C3H6_prop],x=[0.2,0.3,0.5])
#print(f'p={dms_mixed_H2S[0]},c={dms_mixed_H2S[1]},c_err={dms_mixed_H2S[2]}')

# CASE 7: 3-gas mixed Virial & PR DMS : conditional pass on fugacity
dms_mixed_H2S = dms_mixed(press, [H2S_prop,C3H8_prop,C3H6_prop],x=[0.2,0.3,0.5], eos=['pr','virial','virial'],T=[308, 308, 308],gases=['H2S','C3H8', 'C3H6'])
#print(f'p={dms_mixed_H2S[0]},c={dms_mixed_H2S[1]},c_err={dms_mixed_H2S[2]}')

# CASE 8: default selectivity : PASS
dms_mixed_H2S = dms_mixed(press, [H2S_prop,C3H8_prop],x=[0.4,0.6], eos = ['pr','virial'],T=[308,308],gases=['H2S','C3H8'])
dms_mixed_C3H8 = dms_mixed(press, [C3H8_prop,H2S_prop],x=[0.6,0.4], eos = ['virial','pr'],T=[308,308],gases=['C3H8','H2S'])
#print(f'H2S: p={dms_mixed_H2S[0]},c={dms_mixed_H2S[1]},c_err={dms_mixed_H2S[2]}')
#print(f'C3H8: p={dms_mixed_C3H8[0]},c={dms_mixed_C3H8[1]},c_err={dms_mixed_C3H8[2]}')
H2S_C3H8_selec = selectivity(dms_mixed_H2S,dms_mixed_C3H8)
#print(H2S_C3H8_selec)

# CASE 9: custom selectivity : TBD
dms_mixed_H2S = dms_mixed(press, [H2S_prop,C3H8_prop,C3H6_prop],x=[0.2,0.3,0.5], eos=['pr','virial','virial'],T=[308, 308, 308],gases=['H2S','C3H8', 'C3H6'])
dms_mixed_C3H8 = dms_mixed(press, [C3H8_prop,H2S_prop,C3H6_prop],x=[0.3,0.2,0.5], eos=['virial','pr','virial'],T=[308, 308, 308],gases=['C3H8','H2S', 'C3H6'])
dms_mixed_C3H6 = dms_mixed(press, [C3H6_prop,H2S_prop,C3H6_prop],x=[0.5,0.2,0.3], eos=['virial','pr','virial'],T=[308, 308, 308],gases=['C3H6','H2S', 'C3H8'])

def custom_selec_func(x, y, z):
    return x/(y+z)

H2S_C3H8_C3H6_selec = selectivity(dms_mixed_H2S,dms_mixed_C3H8,dms_mixed_C3H6,custom=custom_selec_func)
