/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "incompleteGammaFunction.H"

// * * * * * * * * * * * * * * * * Constructors* * * * * * * * * * * * * * * //

Foam::incompleteGammaFunction::incompleteGammaFunction()
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::incompleteGammaFunction::~incompleteGammaFunction()
{}

// * * * * * * * * * * * * * * * * Member Functions* * * * * * * * * * * * * //

Foam::scalar Foam::incompleteGammaFunction::gammP(const scalar a, const Foam::scalar x) 
{
	if (x < 0.0 || a <= 0.0) throw("bad args in gammP");
	if (x == 0.0) return 0.0;
	else if (a >= aSwitch) return gammPapprox(a,x,1)*Foam::exp(gammln(a));
	else if (x < a + 1.0) return gSer(a,x)*Foam::exp(gammln(a));
	else return (1.0 - gcf(a,x))*Foam::exp(gammln(a));
}

Foam::scalar Foam::incompleteGammaFunction::gammQ(const Foam::scalar a, const Foam::scalar x) 
{
	if (x < 0.0 || a <= 0.0) throw("bad args in gammQ");
	if (x == 0.0) return 1.0;
	else if (a >= aSwitch) return gammPapprox(a,x,0)*Foam::exp(gammln(a));
	else if (x < a + 1.0) return (1.0 - gSer(a,x))*Foam::exp(gammln(a));
	else return gcf(a,x)*Foam::exp(gammln(a));
}

Foam::scalar Foam::incompleteGammaFunction::gSer(const Foam::scalar a, const Foam::scalar x) 
{
	Foam::scalar sum,del,ap;

	gln = gammln(a);
	ap = a;
	del = sum = 1.0/a;
	for(;;){
		++ap;
		del *= x/ap;
		sum += del;
		if (Foam::mag(del) < Foam::mag(sum)*eps){
			return sum*Foam::exp(-x+a*Foam::log(x)-gln);
		}
	}
}

Foam::scalar Foam::incompleteGammaFunction::gcf(const Foam::scalar a, const Foam::scalar x) 
{
	Foam::scalar an,b,c,d,del,h;

	gln = gammln(a);
	b = x + 1.0 - a;
	c = 1.0/fpMin;
	d = 1.0/b;
	h = d;
	
	for (int i=1;;++i) {
		an = -i*(i-a);
		b += 2.0;
		d = an*d + b;
		if (Foam::mag(d) < fpMin) d=fpMin;
		c = b + an/c;
		if (Foam::mag(c) < fpMin) c=fpMin;
		d = 1.0/d;
		del = d*c;
		h *= del;
		if (Foam::mag(del-1.0) <= eps) break;
	}
	return Foam::exp(-x+a*Foam::log(x)-gln)*h;
}

Foam::scalar Foam::incompleteGammaFunction::gammPapprox(Foam::scalar a, Foam::scalar x, Foam::scalar psig)
{
	Foam::scalar xu,t,sum,ans;
	Foam::scalar a1 = a-1.0; 
	Foam::scalar lna1 = Foam::log(a1); 
	Foam::scalar sqrta1 = Foam::sqrt(a1);
	const scalar y[18] = {0.0021695375159141994,
	0.011413521097787704,0.027972308950302116,0.051727015600492421,
	0.082502225484340941, 0.12007019910960293,0.16415283300752470,
	0.21442376986779355, 0.27051082840644336, 0.33199876341447887,
	0.39843234186401943, 0.46931971407375483, 0.54413605556657973,
	0.62232745288031077, 0.70331500465597174, 0.78649910768313447,
	0.87126389619061517, 0.95698180152629142};
	const scalar w[18] = {0.0055657196642445571,
	0.012915947284065419,0.020181515297735382,0.027298621498568734,
	0.034213810770299537,0.040875750923643261,0.047235083490265582,
	0.053244713977759692,0.058860144245324798,0.064039797355015485,
	0.068745323835736408,0.072941885005653087,0.076598410645870640,
	0.079687828912071670,0.082187266704339706,0.084078218979661945,
	0.085346685739338721,0.085983275670394821};
	
	gln = gammln(a);

	if (x > a1) xu = Foam::max(a1 + 11.5*sqrta1, x + 6.0*sqrta1);
	else xu = Foam::max(0.0, Foam::min(a1 - 7.5*sqrta1, x - 5.0*sqrta1));
	sum = 0;
	for (int j=0;j<Foam::incompleteGammaFunction::ngau;++j){
		t = x + (xu-x)*y[j];
		sum += w[j]*Foam::exp(-(t-a1)+a1*(Foam::log(t)-lna1));
	}
	ans = sum*(xu-x)*Foam::exp(a1*(lna1-1.0)-gln);
	return (psig?(x>a1? 1.0-ans:-ans):(x>a1? ans:1.0+ans));
}

Foam::scalar Foam::incompleteGammaFunction::gammln(const scalar xx)
{
	Foam::scalar x,tmp,y,ser;
	static const Foam::scalar cof[14] = {57.1562356658629235,-59.5979603554754912,
								   14.1360979747417471,-0.491913816097620199,.339946499848118887e-4,
								     .465236289270485756e-4,-.983744753048795646e-4,.158088703224912494e-3,
								    -.210264441724104883e-3,.217439618115212643e-3,-.164318106536763890e-3,
									 .844182239838527433e-4,-.261908384015814087e-4,.368991826595316234e-5};

	if (xx <= 0) throw("bad arg in gammln");
	y=x=xx;
	tmp = x+5.24218750000000000;
	tmp = (x+0.5)*Foam::log(tmp)-tmp;
	ser = 0.999999999999997092;
	for (int j=0;j<14;++j) ser += cof[j]/++y;
	return tmp+Foam::log(2.5066282746310005*ser/x);
}

// ************************************************************************* //
