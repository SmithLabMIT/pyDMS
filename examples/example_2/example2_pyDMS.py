import pyDMS.dms as dms

co2 = dms.load_gas_class('example1_pyDMS.pkl')

print(f'\nThe value of b at {co2.temp[0]} K is {co2.b[0]:.4f} atm^-1\n')
print('These results come from the stored arrays of:')
print(f'    temperature: {co2.temp}')
print(f'    b: {co2.b}\n')
