NEURON {
	SUFFIX bk_gp
	USEION k READ ek WRITE ik
	USEION ca READ cai
	RANGE gbar, ik
	GLOBAL minf, mtau, hinf, htau, zinf, ztau
	GLOBAL m_vh, m_k, mtau_y0, mtau_vh1, mtau_vh2, mtau_k1, mtau_k2
	GLOBAL z_coef, ztau
	GLOBAL h_y0, h_vh, h_k, htau_y0, htau_vh1, htau_vh2, htau_k1, htau_k2
	GLOBAL Cq10
}

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(mM) = (milli/liter)
	(S) = (siemens)
}

PARAMETER {
	gbar = 1		(S/cm2)

	m_vh = -28.9		(mV)
	m_k = 6.2		(mV)
	mtau_y0 = .000505	(ms)
	mtau_vh1 = -33.3	(mV)
	mtau_k1 = -10		(mV)
	mtau_vh2 = 86.4		(mV)
	mtau_k2 = 10.1		(mV)

	z_coef = .001		(mM)
	ztau = 1		(ms)

	h_y0 = .085
	h_vh = -32		(mV)
	h_k = 5.8		(mV)
	htau_y0 = .0019		(ms)
	htau_vh1 = -54.2	(mV)
	htau_k1 = -12.9		(mV)
	htau_vh2 = 48.5		(mV)
	htau_k2 = 5.2		(mV)

	cai			(mM)
	celsius	(degC)
	
	Cq10 = 3
}

ASSIGNED {
	g	(S/cm2)
	minf
	mtau	(ms)
	hinf
	htau	(ms)
	zinf
	v	(mV)
	ek	(mV)
	ik	(mA/cm2)
}

STATE {
	m   FROM 0 TO 1
	z   FROM 0 TO 1
	h   FROM 0 TO 1
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar * m * m * m * z * z * h
	ik =  g * (v - ek)
}

DERIVATIVE states {
        rates(v)
        m' = (minf - m) / mtau
        h' = (hinf - h) / htau
        z' = (zinf - z) / ztau
}

PROCEDURE rates(Vm (mV)) {
	LOCAL v, q10
	q10 = Cq10^((celsius - 22 (degC))/10 (degC) )
	v = Vm + 5
	minf = 1 / (1 + exp(-(v - (m_vh)) / m_k))
	mtau = q10*(mtau_y0 + 1 (ms) /(exp((v+ mtau_vh1)/mtau_k1) + exp((v+mtau_vh2)/mtau_k2)))
	zinf = 1/(1 + z_coef / cai)
	hinf = h_y0 + (1-h_y0) / (1+exp((v - h_vh)/h_k))
	htau = q10*(htau_y0 + 1 (ms) /(exp((v + htau_vh1)/htau_k1)+exp((v+htau_vh2)/htau_k2)))
}

INITIAL {
	rates(v)
        m = minf
        z = zinf
        h = hinf
}
