import dadi
def generated_model((n1,n2), pts):
	xx = yy = dadi.Numerics.default_grid(pts)
	phi = dadi.PhiManip.phi_1D(xx)
	m12, m21 = 0, 0
	phi = dadi.Integration.one_pop(phi, xx, 19998.9059512, 13158.5907998)
	phi = dadi.Integration.one_pop(phi, xx, 1.09404879534, 23710.8)
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
	sfs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,yy))
	return sfs
