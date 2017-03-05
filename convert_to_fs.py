import dadi

dd = dadi.Misc.make_data_dict('fin_result_dadi.txt')
fs = dadi.Spectrum.from_data_dict(dd, ['Tan','Nam'], [6,8])
dadi.Spectrum.to_file(fs, 'data.fs')