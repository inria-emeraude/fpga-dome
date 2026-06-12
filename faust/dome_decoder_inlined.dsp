
import("stdfaust.lib");

//////////////// GRIDS /////////////////

cart2spher(x, y, z) = (r, t, p)
                    with{
                    rtemp = sqrt(x^2 + y^2 + z^2);
                    r = rtemp, 1 : select2(rtemp == 0); // avoids r=0
                    t = atan2(y, x);
                    p = asin(z/r);
                    };

/////////// RADIAL ////////////////////

c = 340;
w(r) = 0.5 * c / (r * ma.SR);

secF(1, 1)   =   (1, 1, 0);

secF(2, 1)   =   (1, 3, 3);

secF(3, 1)   =   (1, 2.3221853546260855929, 0);
secF(3, 2)   =   (1, 3.6778146453739144071, 6.4594326934833653473);

secF(4, 1)   =   (1, 4.2075787943592556632, 11.4878004768711997988);
secF(4, 2)   =   (1, 5.7924212056407443368, 9.1401308902779310256);

secF(5, 1)   =   (1, 3.6467385953296432597, 0);
secF(5, 2)   =   (1, 6.7039127983070662860, 14.272480513279948265);
secF(5, 3)   =   (1, 4.6493486063632904542, 18.156315313452237137);

secF(6, 1)   =   (1, 8.4967187917267278899, 18.801130589570517411);
secF(6, 2)   =   (1, 7.4714167126516293359, 20.852823177396347991);
secF(6, 3)   =   (1, 5.0318644956216427742, 26.514025344068052456);

secF(7, 1)   =   (1, 4.9717868585279356779, 0);
secF(7, 2)   =   (1, 9.5165810563092578905, 25.666444752769034175);
secF(7, 3)   =   (1, 8.1402783272762749434, 28.936546093263966238);
secF(7, 4)   =   (1, 5.3713537578865314883, 36.596785156877450848);

secF(8, 1)   =   (1, 11.1757720865261703980, 31.977225258279201354);
secF(8, 2)   =   (1, 10.4096815812737638365, 33.934740085181713765);
secF(8, 3)   =   (1, 8.7365784344048048141, 38.569253275096191935);
secF(8, 4)   =   (1, 5.6779678977952609514, 48.432018652637095880);

secF(9, 1)   =   (1, 6.2970191817149685378, 0);
secF(9, 2)   =   (1, 12.2587358085485455756, 40.589267909914637799);
secF(9, 3)   =   (1, 11.2088436390155628324, 43.646645753129244892);
secF(9, 4)   =   (1, 9.2768797743607805933, 49.788502657376288447);
secF(9, 5)   =   (1, 5.9585215963601424609, 62.041437621985133043);

secF(10, 1)  =   (1, 13.8440898108544922308, 48.667548564148698918);
secF(10, 2)  =   (1, 13.2305819309537405179, 50.582361562872006750);
secF(10, 3)  =   (1, 11.9350566571755716807, 54.839156202307484983);
secF(10, 4)  =   (1, 9.7724391337179991598, 62.625585912537518586);
secF(10, 5)  =   (1, 6.2178324672981964107, 77.442700531277433593);

tf2sum(b0, b1, b2, a1, a2) = sub~sum1(a1, a2): sum2(b0, b1, b2)
	with {
		sum1(k1, k2, x) 	= x:(+~_<:((_':+~_),*(k1)):*(k2),_:+);
		sum2(k0, k1, k2, x)  =x<:*(k0),+~_,_:_,(-<:*(k1),(_':+~_)*(k2):+):+;
		sub(x, y)		= y - x;
	};

nfc(0,r) = _ * r;
nfc(l,r) = seq(i, ceil(l / 2), secNFC(l, i+1, r)) * r;

secNFC(l,numsec,r) = _*1/g2(r):tf2sum(1,0,0,d21(r),d22(r))
    with {
        bp(r)   =   secF(l,numsec):(_*1,_*w(r),_*w(r)^2); // b' Coefficients
        g2(r)   =   bp(r):(+,_):+; // g2 Coefficients
        d21(r)  =   bp(r):(!,_*2,_*4):+/g2(r); // d21 Coefficients
        d22(r)  =   bp(r):(!,!,_*4)/g2(r); // d22 Coefficients
    };

/////////// YLM ////////////////////////

wre(lmax, l1, l) = legendre(l,cos(137.9 * ma.PI / 180 / (l1 + 1.51))) / sum(ll, lmax + 1, (2 * ll + 1) * legendre(ll, cos(137.9 * ma.PI / 180 / (l1 + 1.51))) * (ll <= l1)) * (l <= l1);

legendre(l,x) = case{
                (0) => 1;
                (1) => x;
                (l) => ((2*l-1)*x*legendre(l-1,x) - (l-1)*legendre(l-2,x))/l;
                }(l);

ylm(l,m,t,p) = n3d(l,m)*alegendre(l,abs(m),sin(p))*
            case{
            //(1) => sin(abs(m)*t);
            (1) => chebyshev2(abs(m)-1,cos(t))*sin(t); // [8];
            //(0) => cos(abs(m)*t);
            (0) => ma.chebychev(abs(m),cos(t)); // [8]
            }(m<0);

alegendre(l,m,x) =  case{
                (1,0,1) => 1; // special case to avoid evaluate factorial2(-1)
                (1,0,0) => factorial2(2*l-1)*(1-x^2)^(l/2); // (1)^m not included here
                (0,1,0) => x*(2*l-1)*alegendre(l-1,l-1,x);
                (0,0,0) => 1/(l-m)*((2*l-1)*x*alegendre(l-1,m,x)-(l-1+m)*alegendre(l-2,m,x));
                }(m==l,m==(l-1),l==0);

factorial(m) = ma.gamma(m+1);

factorial2(m) = 2^(m/2 + 1/4*(1 - cos(m*ma.PI)))*ma.PI^(1/4*(-1 + cos(m*ma.PI)))*ma.gamma(1+m/2);

n3d(l,m) =  sqrt((2*l+1)*factorial(l-abs(m))/factorial(l+abs(m)))*
            case{
            (0) => 1;
            (m) => sqrt(2);
            }(m);

chebyshev2(m,x) = case{
                (0) => 1;
                (1) => 2*x;
                (m) => 2*x*chebyshev2(m-1,x) - chebyshev2(m-2,x);
                }(m);

sup(c) = R(c) with {
 R((c,cl)) = max(R(c),R(cl));
 R(c)      = c;
};

// COMPILATION PARAMETERS
L	=	3; // Maximum degree $L$.
N   =   17; // Loudspeaker number
nfcon = 1; // Include NFC (1=included, 0=not included)

// Loudspeakers Cartesian coordinates $(x, y, z)$ in meters.
speaker(0) = (0.71, 0, 0);
speaker(1) = (0.502, 0.502, 0);
speaker(2) = (0, 0.71, 0);
speaker(3) = (0.502, -0.502, 0);
speaker(4) = (-0.71, 0, 0);
speaker(5) = (-0.502, -0.502, 0);
speaker(6) = (0, -0.71, 0);
speaker(7) = (0.502, -0.502, 0);
speaker(8) = (0.4, 0, 0.62);
speaker(9) = (0.283, 0.283, 0.62);
speaker(10) = (0, 0.4, 0.62);
speaker(11) = (0.283, -0.283, 0.62);
speaker(12) = (0, -0.4, 0.62);
speaker(13) = (-0.283, -0.283, 0.62);
speaker(14) = (0, -0.4, 0.62);
speaker(15) = (0.283, -0.283, 0.62);
speaker(16) = (0, 0, 0.79);

// DO NOT EDIT BELOW HERE

// Inputs/Outputs
ins	=	(L+1)^2;
outs	=	N;

rs = par(i, N, ba.take(1, cart2spher(ba.take(1, speaker(i)), ba.take(2, speaker(i)), ba.take(3, speaker(i)))));
rmax = sup(rs);

// User Interface
volout	=	hslider("Outputs Gain[unit:dB][style:knob][osc:/levelout -70 6]", 0, -70, 6, 0.1) : ba.db2linear;// : si.smoo;

row(i) = case{
                (0) => par(l, L + 1, par(m, 2 * l + 1, _ * ylm(l, m - l, t, p)) :> nfc(l,r) * wre(L, L, l) ) :>_ / rmax : de.delay(rmax/c, (rmax-r)/c); // normalize by rmax
                (1) => par(l, L + 1, par(m, 2 * l + 1, _ * ylm(l, m - l, t, p)) :> _* wre(L, L, l)) :>_;
                }(nfcon==0)
                with {
                    r = ba.take(1, cart2spher(ba.take(1, speaker(i)), ba.take(2, speaker(i)), ba.take(3, speaker(i))));
                    t = ba.take(2, cart2spher(ba.take(1, speaker(i)), ba.take(2, speaker(i)), ba.take(3, speaker(i))));
                    p = ba.take(3, cart2spher(ba.take(1, speaker(i)), ba.take(2, speaker(i)), ba.take(3, speaker(i))));
                    };

matrix	=	si.bus(ins)<:par(i, N, row(i) * (volout));

process = matrix;
