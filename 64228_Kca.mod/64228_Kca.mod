:Calcium Activated Potassium Channels

NEURON 
{
	SUFFIX Kca
	:USEION Ca READ Cai VALENCE 2
	:USEION Kca WRITE iKca VALENCE 1
	USEION ca READ cai
        USEION k READ ek WRITE ik
        RANGE infmKcaV,taumKcaV,gKcabar,Cahalf
}

UNITS
{
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mS) = (millimho)
	(mol)= (1)
	(M)  = (mol/liter)
	(uM) = (micro M)
}

PARAMETER
{
       
 
       :Ca-dependent K current
       :eKca=-80 (mV)
       gKcabar = 5 (mS/cm2)
       :Cai
	Cahalf=1 (uM)	 
 

}

STATE
{

	:mKcaV
	mKcaCa
	
}

ASSIGNED
{
        ek (mV)
        cai (mM)
	ik (mA/cm^2)
	v (mV)
           :Ca-dependent potassium channel, Kca
	:infmKcaV
	:taumKcaV  (ms)
	gKca (mho/cm2)

}

INITIAL
{      LOCAL Cas
	:rate(v)
	Cas=cai*1000 :uM

	:mKcaV= infmKcaV
        :Cas=cai
        mKcaCa=1/(1+(Cahalf/Cas)^4)
}




BREAKPOINT
{
	LOCAL Cas
	:SOLVE states METHOD cnexp
	Cas=cai*1000 :uM
	mKcaCa=1/(1+(Cahalf/Cas)^4 )
	:gKca=(0.001)*gKcabar*(mKcaV^2)*mKcaCa
        gKca=(0.001)*gKcabar*mKcaCa^4
	ik=gKca*(v-ek) 
	: the current is in the unit of mA/cm2
	
}

