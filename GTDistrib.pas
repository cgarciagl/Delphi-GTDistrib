unit GTDistrib;

interface

uses  Windows, Messages, SysUtils, Classes, Graphics, Controls, Forms, Dialogs, math;

const
  eps                   = 5E-12;
  maxparametercount     = 3;
  maxvalue              = maxlongint;
  id_error              = 1300;
  err_val_le	          = 1301;
  err_val_l	            = 1302;
  err_val_ge	          = 1303;
  err_val_int           = 1304;
  err_val2_le	          = 1305;
  err_val2_l	          = 1306;
  err_val2_ge	          = 1307;
  err_val2_int          = 1308;
  err_val3_le	          = 1309;
  err_val3_l	          = 1310;
  err_val3_ge	          = 1311;
  err_val3_int          = 1312;
  err_val1_l_val2       = 1313;
  err_val1_l_val3       = 1314;
  err_val_g             = 1315;
  err_overflow          = 1361;
  err_expectation       = 1362;
  err_variance          = 1363;

var seed1,seed2:integer;

type

EArgumentError = class(Exception);

range = record
  lower, upper:extended;
end;

parameter = record
  value: extended;
  integertype: boolean;
  bound: range;
end;

parameterarray = array[1..maxparametercount] of parameter;

tDistribution = class
private
  Fparametercount: 0..maxparametercount;
  Fparameters: parameterarray;
  support: range;
  function parameter_ok(i:integer):boolean;
public
  constructor Create;
  destructor free;
  function Setparameters(par1,par2,par3:extended):boolean; virtual;
  function in_support(x:extended):boolean; virtual;
  function density(x:extended):extended; virtual;
  function probability(x:extended): extended; virtual;
  function cdf(x:extended):extended; virtual;
  function prob_between(a,b:extended):extended; virtual;
  function quantile(var q:extended):extended; virtual; abstract;
  function expectation_exists:boolean;virtual;
  function variance_exists:boolean;virtual;
  function Mean:extended; virtual; abstract;
  function Variance:extended; virtual; abstract;
  function SD:extended; virtual;
  function random_value:extended; virtual;
end;

tContinuous_distribution = class(tDistribution)
public
  function quantile(var q:extended):extended; override;
end;

tNormal_distribution = class(tContinuous_distribution)
private
  mu, sigma2, sigma, const_part, random_const: extended;
public
  function Setparameters(par1,par2,par3:extended):boolean; override;
  constructor create(newmu, newsigma2:extended);
  function density(x:extended):extended; override;
  function cdf(x:extended):extended; override;
  function random_value:extended; override;
  function Mean:extended; override;
  function Variance:extended; override;
end;

tBeta_distribution = class(tContinuous_distribution)
private
  alpha, beta, const_part, alpha_plus_beta, random_const1, random_const2: extended;
public
  function Setparameters(par1,par2,par3:extended):boolean; override;
  constructor create(newalpha, newbeta:extended);
  function in_support(x:extended):boolean; override;
  function density(x:extended):extended; override;
  function cdf(x:extended):extended; override;
  function random_value:extended; override;
  function Mean:extended; override;
  function Variance:extended; override;
end;

tGamma_distribution = class(tContinuous_distribution)
private
  alpha, lambda, const_part, random_const1, random_const2: extended;
public
  function Setparameters(par1,par2,par3:extended):boolean; override;
  constructor create(newalpha, newlambda:extended);
  function in_support(x:extended):boolean; override;
  function density(x:extended):extended; override;
  function cdf(x:extended):extended; override;
  function random_value:extended; override;
  function Mean:extended; override;
  function Variance:extended; override;
end;

tChi_square_distribution = class(tGamma_distribution)
private
  df: extended;
public
  function Setparameters(par1,par2,par3:extended):boolean; override;
  constructor create(newdf:extended);
end;

tExponential_distribution = class(tGamma_distribution)
private
  theta: extended;
public
  function Setparameters(par1,par2,par3:extended):boolean; override;
  constructor create(newtheta:extended);
  function quantile(var q:extended):extended; override;
end;

tF_distribution = class(tContinuous_distribution)
private
  df1, df2, const_part: extended;
  cbeta: tBeta_distribution;
  cp: tChi_square_distribution;
public
  function Setparameters(par1,par2,par3:extended):boolean; override;
  constructor create(newdf1,newdf2:extended);
  destructor free;
  function in_support(x:extended):boolean; override;
  function density(x:extended):extended; override;
  function cdf(x:extended):extended; override;
  function random_value:extended; override;
  function expectation_exists:boolean; override;
  function variance_exists:boolean; override;
  function Mean:extended; override;
  function Variance:extended; override;
end;

tT_distribution = class(tContinuous_distribution)
private
  df, const_part, random_const: extended;
  N: tNormal_distribution;
  G: tGamma_distribution;
public
  function Setparameters(par1,par2,par3:extended):boolean; override;
  constructor create(newdf:extended);
  destructor free;
  function density(x:extended):extended; override;
  function cdf(x:extended):extended; override;
  function random_value:extended; override;
  function expectation_exists:boolean; override;
  function variance_exists:boolean; override;
  function Mean:extended; override;
  function Variance:extended; override;
end;

tCauchy_distribution = class(tContinuous_distribution)
private
  alpha, beta:extended;
public
  function Setparameters(par1,par2,par3:extended):boolean; override;
  constructor create(newalpha, newbeta:extended);
  function density(x:extended):extended; override;
  function cdf(x:extended):extended; override;
  function random_value:extended; override;
  function expectation_exists:boolean;override;
  function variance_exists:boolean;override;
end;

tDouble_exponential_distribution = class(tContinuous_distribution)
private
  alpha, beta:extended;
public
  function Setparameters(par1,par2,par3:extended):boolean; override;
  constructor create(newalpha, newbeta:extended);
  function density(x:extended):extended; override;
  function cdf(x:extended):extended; override;
  function random_value:extended; override;
  function Mean:extended; override;
  function Variance:extended; override;
end;

tGumbel_distribution = class(tContinuous_distribution)
private
  alpha, beta:extended;
public
  function Setparameters(par1,par2,par3:extended):boolean; override;
  constructor create(newalpha, newbeta:extended);
  function density(x:extended):extended; override;
  function cdf(x:extended):extended; override;
  function random_value:extended; override;
  function Mean:extended; override;
  function Variance:extended; override;
end;

tLogistic_distribution = class(tContinuous_distribution)
private
  alpha, beta:extended;
public
  function Setparameters(par1,par2,par3:extended):boolean; override;
  constructor create(newalpha, newbeta:extended);
  function density(x:extended):extended; override;
  function cdf(x:extended):extended; override;
  function random_value:extended; override;
  function Mean:extended; override;
  function Variance:extended; override;
end;

tLognormal_distribution = class(tContinuous_distribution)
private
  mu, sigma2, sigma, const_part: extended;
  norm: tNormal_distribution;
public
  function Setparameters(par1,par2,par3:extended):boolean; override;
  constructor create(newmu, newsigma2:extended);
  destructor free;
  function in_support(x:extended):boolean; override;
  function density(x:extended):extended; override;
  function cdf(x:extended):extended; override;
  function random_value:extended; override;
  function Mean:extended; override;
  function Variance:extended; override;
end;

tPareto_distribution = class(tContinuous_distribution)
private
  alpha, theta:extended;
public
  function Setparameters(par1,par2,par3:extended):boolean; override;
  constructor create(newalpha, newtheta:extended);
  function density(x:extended):extended; override;
  function cdf(x:extended):extended; override;
  function random_value:extended; override;
  function expectation_exists:boolean; override;
  function variance_exists:boolean; override;
  function Mean:extended; override;
  function Variance:extended; override;
end;

tUniform_distribution = class(tContinuous_distribution)
private
  alpha, beta:extended;
public
  function Setparameters(par1,par2,par3:extended):boolean; override;
  constructor create(newalpha, newbeta:extended);
  function density(x:extended):extended; override;
  function cdf(x:extended):extended; override;
  function random_value:extended; override;
  function Mean:extended; override;
  function Variance:extended; override;
end;

tWeibull_distribution = class(tContinuous_distribution)
private
  a, b, b1:extended;
public
  function Setparameters(par1,par2,par3:extended):boolean; override;
  constructor create(newa, newb:extended);
  function in_support(x:extended):boolean; override;
  function density(x:extended):extended; override;
  function cdf(x:extended):extended; override;
  function random_value:extended; override;
  function Mean:extended; override;
  function Variance:extended; override;
end;

tNon_central_chi_square_distribution = class(tContinuous_distribution)
private
  df, lambda:extended;
  central_chi2: tChi_square_distribution;
  c_chi2: tChi_square_distribution;
  stand_norm: tNormal_distribution;
public
  function Setparameters(par1,par2,par3:extended):boolean; override;
  constructor create(newdf,newlambda:extended);
  destructor free;
  function density(x:extended):extended; override;
  function cdf(x:extended):extended; override;
  function Mean:extended; override;
  function Variance:extended; override;
end;

tNon_central_Beta_distribution = class(tContinuous_distribution)
private
  alpha, beta,lambda:extended;
public
  function Setparameters(par1,par2,par3:extended):boolean; override;
  constructor create(newalpha,newbeta,newlambda:extended);
  function in_support(x:extended):boolean; override;
  function density(x:extended):extended; override;
  function cdf(x:extended):extended; override;
  function Mean:extended; override;
  function Variance:extended; override;
end;

tNon_central_F_distribution = class(tContinuous_distribution)
private
  df1, df2, lambda:extended;
  ncBeta: tNon_central_Beta_distribution;
public
  function Setparameters(par1,par2,par3:extended):boolean; override;
  constructor create(newdf1,newdf2,newlambda:extended);
  destructor free;
  function in_support(x:extended):boolean; override;
  function density(x:extended):extended; override;
  function cdf(x:extended):extended; override;
  function expectation_exists:boolean;override;
  function variance_exists:boolean;override;
  function Mean:extended; override;
  function Variance:extended; override;
end;

tNon_central_t_distribution = class(tContinuous_distribution)
private
  df, delta: extended;
  central_t: tT_distribution;
public
  function Setparameters(par1,par2,par3:extended):boolean; override;
  constructor create(newdf,newdelta:extended);
  destructor free;
  function density(x:extended):extended; override;
  function cdf(x:extended):extended; override;
  function expectation_exists:boolean;override;
  function variance_exists:boolean;override;
  function Mean:extended; override;
  function Variance:extended; override;
end;

tDiscrete_distribution = class(tDistribution)
private
  mode, maxprob, random_const:extended;
public
  procedure fill_constants; virtual;
  function in_support(x:extended):boolean; override;
  function factor(i:integer):extended; virtual;
         {  ratio probability(i)/probability(i-1) in order to speed up
            computation of probabilities in prob_between and quantile     }
  function prob_between(a,b:extended):extended; override;  { P(a <= X <= b) }
  function quantile(var q:extended):extended; override;
  function random_value:extended; override;
end;

tBinomial_distribution = class(tDiscrete_distribution)
private
  n: integer; p: extended;
public
  function Setparameters(par1,par2,par3:extended):boolean; override;
  constructor create(newn, newp:extended);
  function factor(i:integer):extended; override;
  function probability(x:extended): extended; override;
  function cdf(x:extended):extended; override;
  function Mean:extended; override;
  function Variance:extended; override;
end;

tBernoulli_distribution = class(tBinomial_distribution)
public
  constructor create(newp:extended);
  function quantile(var q:extended):extended; override;
  function random_value:extended; override;
end;

tHypergeometric_distribution = class(tDiscrete_distribution)
private
  popul_N, N1, sample_n, N2: integer;
public
  function Setparameters(par1,par2,par3:extended):boolean; override;
  constructor create(newpopul_N,newN1,newsample_n:extended);
  function probability(x:extended): extended; override;
  function factor(i:integer):extended; override;
  function Mean:extended; override;
  function Variance:extended; override;
end;

tPoisson_distribution = class(tDiscrete_distribution)
private
  lambda: extended;
public
  function Setparameters(par1,par2,par3:extended):boolean; override;
  constructor create(newlambda:extended);
  function probability(x:extended):extended; override;
  function cdf(x:extended):extended; override;
  function Mean:extended; override;
  function Variance:extended; override;
end;

tDiscrete_uniform_distribution = class(tDiscrete_distribution)
private
  N: integer;
public
  function Setparameters(par1,par2,par3:extended):boolean; override;
  constructor create(newN:extended);
  function probability(x:extended):extended; override;
  function random_value:extended; override;
  function Mean:extended; override;
  function Variance:extended; override;
end;

tNegative_binomial_distribution = class(tDiscrete_distribution)
private
  r: integer; p: extended;
public
  function Setparameters(par1,par2,par3:extended):boolean; override;
  constructor create(newr,newp:extended);
  function probability(x:extended):extended; override;
  function factor(i:integer):extended; override;
  function Mean:extended; override;
  function Variance:extended; override;
end;

tGeometric_distribution = class(tNegative_binomial_distribution)
public
  constructor create(newp:extended);
end;

implementation

{-------------------------------------------------------------------}
{  Auxiliary Functions:                                             }
{-------------------------------------------------------------------}

function approx(x, y:extended):Boolean;
Begin
  result:= abs(x-y) < eps;
End;

Function nearint(x: extended): Boolean;
Begin
  result:= (Abs(x) < maxlongint) and approx(x,round(x));
End;

function ExtToStr(x:extended):string;
begin
  result:=FormatFloat('0.############',x);
end;

function Check(condition:boolean; err_nr:integer):boolean;
begin
  if condition
  then raise EArgumentError.Create(LoadStr(err_nr));
  result:=not condition;
end;

function Checkval(val1:extended; err_nr:integer; val2:extended):boolean;
begin
  result:=false;
  case err_nr of
  err_val_l,err_val2_l,err_val3_l:    result:=val1<val2;
  err_val_le,err_val2_le,err_val3_le: result:=val1<=val2;
  err_val_g:                          result:=val1>val2;
  err_val_ge,err_val2_ge,err_val3_ge: result:=val1>=val2;
  end;
  if result then raise EArgumentError.Create(LoadStr(err_nr)+ExtToStr(val2));
  result:=not result;
end;


Function factorial(n: extended): extended;
Var i:integer;
Begin
  result:= 1;
  for i:=2 to round(n) do result:=result/i;
  if result*maxlongint>1 then result:=round(1/result)
  else result:=1/result;
End;

Function log_choose(n, m: extended): extended;
Var i: integer;
Begin
  result:=0;
  if Checkval(n,err_val_le,0) and Check(not nearint(n),err_val_int) and Check(n>maxlongint,err_overflow)
  and Checkval(m,err_val_l,0) and Check(not nearint(m),err_val_int) and Check(m>n,err_val1_l_val2)
  then
    begin
      If n - m < m Then m := n - m;
      For i := 1 To Round(m) Do result := result + Ln((n - i + 1) / i);
    end;
End;

function Floor(x: Extended): Integer;
begin
  Result := Trunc(x);
  if Frac(x) < 0 then Dec(Result);
end;

Function log_gamma(x:extended):extended; {Numerical Recipes in Pascal}
var y,tmp,ser:extended; j:integer; cof:array[1..6] of extended;
begin
  cof[1]:=76.18009172947146; cof[2]:=-86.50532032941677;
  cof[3]:=24.01409824083091; cof[4]:=-1.231739572450155;
  cof[5]:=0.1208650973866179E-2; cof[6]:=-0.5395239384953E-5;
  y:=x-1; tmp:=y+5.5; tmp:=(y+0.5)*ln(tmp)-tmp; ser:=1;
  for j:=1 to 6 do
    begin
      y:=y+1; ser:=ser+cof[j]/y;
    end;
    result:=tmp+ln(2.5066282746310005*ser);
end;

Function incgamma(x, a: extended):extended; {Numerical Recipes In Pascal,p.181}

function gser(x, a:extended):extended; { Numerical Recipes In Pascal, p.181}
var n:integer; sum,del,ap:extended;
begin
  if approx(0,x) then result:=0
  else
    begin
      ap:=a; sum:=1/a; del:=sum; n:=1;
      repeat
        ap:=ap+1;
        del:=del*x/ap;
        sum:=sum+del;
        inc(n);
      until (abs(del) < abs(sum)*eps) or (n>100);
      result:=sum*exp(-x+a*ln(x)-log_gamma(a));
    end;
end;

function gcf(x, a:extended):extended; { Numerical Recipes In Pascal, p.181}
var n:integer; gold,g,fac,b1,b0,anf,ana,a1,a0:extended;
begin
  g:=0; gold:=1; a0:=1; a1:=x; b0:=0; b1:=1; fac:=1; n:=1;
  repeat
    ana:=n-a;
    a0:=(a1+a0*ana)*fac; b0:=(b1+b0*ana)*fac; anf:=n*fac;
    a1:=x*a0+anf*a1; b1:=x*b0+anf*b1;
    if a1<>0 then
      begin
        gold:=g;
        fac:=1/a1;
        g:=b1*fac;
      end;
    inc(n);
  until abs((g-gold)/g)<eps;
  result:=exp(-x+a*ln(x)-log_gamma(a))*g;
end;

begin  {incgamma}
  result:=1;
  if Checkval(x,err_val_l,0) and Check(x>maxvalue,err_overflow)
  and Checkval(a,err_val_le,0) and Check(a>maxvalue,err_overflow)
  then
  if x < a+1 then result:= gser(x,a)
  else result:= 1 - gcf(x,a);
end;   {incgamma}

Function incbeta(x, a, b: extended): extended;
Var bt:extended;

function betacf(x,a,b:extended):extended;
 { continued fraction (26.5.8) from Abramowitz & Stegun ;
   p.188 from Numerical Recipes in Pascal, Press, Flannery, Teukolsky, Vetterling }
var tem, qap, qam, qab,d, bz, bpp,bp,bm,az,app,  am,aold,ap:extended;
    m:integer;
begin
  am:=1; bm:=1; az:=1;
  qab:=a+b; qap:=a+1; {these q's will be used in factors which occur in the coeffs}
  qam:=a-1; bz:=1-qab*x/qap; m:=1;
  repeat
    tem:=2*m;
    d:=m*(b-m)*x/((qam+tem)*(a+tem));
    ap:=az+d*am; bp:=bz+d*bm;
    d:=-(a+m)*(qab+m)*x/((a+tem)*(qap+tem));
    app:=ap+d*az; bpp:=bp+d*bz;
    aold:=az;
    am:=ap/bpp; bm:=bp/bpp; az:=app/bpp;
    bz:=1; inc(m);
  until abs(az-aold) < eps*abs(az);
  result:=az;
end;

Begin   {incbeta}
try
  If x<=0 then result:=0
  else if x>=1 then result:=1
  else
    begin
      bt:=exp(a*ln(x)+b*ln(1-x)+log_gamma(a+b)-log_gamma(a)-log_gamma(b));
      if x < (a+1)/(a+b+2)
      then result:= bt*betacf(x,a,b)/a
      else result:=1.0 - bt*betacf(1-x,b,a)/b;
    end;
except result:=1;
end;
end;    {incbeta}

Function randomdraw:extended;
{From: L'Ecuyer,P: Efficient and Portable Combined Random Number Generators,
 Comm of the ACM, June 1988, Vol31, Nr 6}
const m1=2147483563; a1=40014; m2=2147483399; a2=40692; q1=53668; q2=52774;
var Z,k:integer;
begin
  k:=seed1 div q1;
  seed1:=a1*(seed1-k*q1) - k*12211;
  if seed1<0 then seed1:=seed1+m1;
  k:=seed2 div q2;
  seed2:=a2*(seed2-k*q2) - k*3791;
  if seed2<0 then seed2:=seed2+m2;
  Z:=seed1-seed2;
  if Z<1 then Z:=Z+m1-1;
  result:=Z*4.656613E-10;
end;

{-------------------------------------------------------------------}
{  tDistribution                                                    }
{-------------------------------------------------------------------}


constructor tDistribution.create;
begin
  Fparametercount:=2;
  support.lower := -maxvalue; support.upper := maxvalue;
end;

destructor tDistribution.free;
begin
  inherited destroy;
end;

function tDistribution.Setparameters(par1,par2,par3:extended):boolean;
var i:1..maxparametercount+1;
begin
  Fparameters[1].value:=par1; Fparameters[2].value:=par2; Fparameters[3].value:=par3;
  i:=1;
  while (i<=Fparametercount) and parameter_ok(i) do inc(i);
  result:=i>Fparametercount;
  support.lower:=-maxvalue; support.upper:=maxvalue;
end;

function tDistribution.parameter_ok(i:integer):boolean;
begin
  with Fparameters[i] do
    begin
      if integertype
      then result:= Checkval(value,{err_val_l}id_error,bound.lower)
        and Checkval(value,{err_val_g}id_error+4*i-1,bound.upper)
        and Check(not nearint(value),{err_val_int}id_error+4*i)
      else result:= Checkval(value,{err_val_le}id_error+4*i-3,bound.lower)
        and Checkval(value,err_val_ge{id_fout+4*i-1},bound.upper);
    end;
end;

function tDistribution.SD:extended;
begin
  Check(not variance_exists,err_variance);
  result:=sqrt(Variance);
end;

function tDistribution.in_support(x:extended):boolean;
begin
 result:=(x>=support.lower) and (x<=support.upper);
end;

function tDistribution.density(x:extended):extended;
begin
  result:=0;
end;

function tDistribution.probability(x:extended): extended;
begin
  result:=0;
end;

function tDistribution.cdf(x:extended):extended;
begin
  if x<support.lower then result:=0
  else
    if x>=support.upper then result:=1
    else result:=prob_between(support.lower,x);
end;

function tDistribution.prob_between(a,b:extended):extended;
var prob:extended;
begin
  if a<support.lower then a:=support.lower;
  if b>support.upper then b:=support.upper;
  if b >= maxvalue then prob := 1 - cdf(a) + probability(a)
  else if a <= -maxvalue then prob := cdf(b)
  else prob := cdf(b) - cdf(a) + probability(a);
  if prob < 0 then result:= 0
  else if prob > 1 then result:= 1
  else result:= prob;
end;

function tDistribution.random_value:extended;
var x:extended;
begin
  x:=randomdraw; result:=quantile(x);
end;

function tDistribution.expectation_exists:boolean;
begin
  result:=true;
end;

function tDistribution.variance_exists:boolean;
begin
  result:=true;
end;

{-------------------------------------------------------------------}
{ tContinuous_Distribution                                          }
{-------------------------------------------------------------------}

function tContinuous_Distribution.quantile(var q:extended):extended;
var pdf,dx,G,xhi,xlo,x,xold:extended;
begin
  if (not Checkval(q,err_val_l,0) or approx(0,q)) then quantile:=support.lower
  else if (not Checkval(q,err_val_g,1) or approx(1,q)) then quantile:=support.upper
  else
  begin
    xlo:=support.lower; xhi:=support.upper;
    if expectation_exists
    then x:=Mean
    else x:=(xlo+xhi)/2;

    if variance_exists
    then dx:=2*SD
    else dx:=8;

    if x+dx<xhi then xhi:=x+dx;
    if x-dx>xlo then xlo:=x-dx;

    while cdf(xhi)<q do
      begin
        xlo:=xhi; xhi:=x+8*(xhi-x);
        if xhi>support.upper then xhi:=support.upper;
      end;
    if approx(xhi-xlo,2*dx)        { the `while' body was not executed }
    then
    while cdf(xlo)>q do
      begin
        xhi:=xlo; xlo:=x+8*(xhi-x);
        if xlo<support.lower then xlo:=support.lower;
      end;
    x:=(xlo+xhi)/2;                { xlo<=x<=xhi and cdf(xlo)<q<cdf(xhi) }

    repeat
      G:=cdf(x)-q; pdf:=density(x); xold:=x;
      if G<=0 then xlo:=x; if G>=0 then xhi:=x;
      if (xhi-xlo)*pdf<=abs(2*G)   { if Newton-Raphson would go outside [xlo,xhi] }
      then
        begin
          x:=(xlo+xhi)/2;          { bisection }
        end
      else
        begin
          x:=x-G/pdf;              { Newton-Raphson }
        end;
    until abs(x-xold) < eps;
    result:=x;
  end;
end;

{-------------------------------------------------------------------}
{ tNormal_distribution                                              }
{-------------------------------------------------------------------}

constructor tNormal_distribution.Create(newmu, newsigma2:extended);
begin
  Fparametercount:=2;
  with Fparameters[1] do
    begin
      bound.lower:=-maxvalue; bound.upper:=maxvalue;
      integertype:=false; value:=newmu;
    end;
  with Fparameters[2] do
    begin
      bound.lower:=0; bound.upper:=maxvalue;
      integertype:=false; value:=newsigma2;
    end;
  Setparameters(newmu,newsigma2,0);
end;

function tNormal_distribution.Setparameters(par1,par2,par3:extended):boolean;
begin
  result:=inherited Setparameters(par1,par2,par3);
  mu:=par1; sigma2:=par2; sigma:=sqrt(sigma2);
  const_part:= -(ln(2*pi)+ln(par2))/2;
  random_const:=sqrt(2/exp(1));
end;

function tNormal_distribution.density(x:extended):extended;
begin
  if in_support(x)
  then result:= exp(const_part-sqr(x-mu)/(2*sigma2))
  else result:=0;
end;

function tNormal_distribution.cdf(x:extended):extended;
var z:extended;
begin
  z:=(x-mu)/sigma;
  if in_support(x)
  then
    if z>=0 then result:=0.5+0.5*incgamma(0.5*sqr(z),0.5)
    else result:=0.5-0.5*incgamma(0.5*sqr(z),0.5)
  else result:=inherited cdf(x);
end;

function tNormal_distribution.random_value:extended;
var x,u,v,x2:extended;
begin  {Devroye: Non-Uniform Random Variate Generation, p.197,199}
  repeat
    u:=randomdraw;
    v:=random_const*(2*randomdraw-1);
    x:=v/u; x2:=sqr(x);
  until (x2<=6-8*u+2*sqr(u)) or ((x2<2/u-2*u) and (x2<=-4*ln(u)));
  result:=mu+sigma*x;
end;

function tNormal_distribution.Mean:extended;
begin
  result:=mu;
end;

function tNormal_distribution.Variance:extended;
begin
  result:=sigma2;
end;

{-------------------------------------------------------------------}
{  tBeta_distribution                                               }
{-------------------------------------------------------------------}

constructor tBeta_distribution.create(newalpha, newbeta:extended);
begin
  Fparametercount:=2;
  With FParameters[1] do
    begin
      bound.lower:=0; bound.upper:=maxvalue;
      integertype:=false;
    end;
  With Fparameters[2] do
    begin
      bound.lower:=0; bound.upper:=maxvalue;
      integertype:=false;
    end;
  Setparameters(newalpha, newbeta, 0);
end;

function tBeta_distribution.Setparameters(par1,par2,par3:extended):boolean;
begin
  result:=inherited Setparameters(par1,par2,par3);
  alpha:=par1; beta:=par2; alpha_plus_beta:=alpha+beta;
  const_part:=log_gamma(alpha_plus_beta)-log_gamma(alpha)-log_gamma(beta);
  support.lower := 0; support.upper :=1;
  if ((alpha<=1) or (beta<=1)) then
    begin
      if alpha<beta then random_const1:=alpha else random_const1:=beta;
    end
  else random_const1:=sqrt((2*alpha*beta-alpha_plus_beta)/(alpha_plus_beta-2));
  random_const2:=alpha+random_const1;
end;

function tBeta_distribution.in_support(x:extended):boolean;
begin
  if (alpha<1) and (beta<1) then result:=(x>0) and (x<1)
  else if alpha<1 then result:=(x>0) and (x<=1)
  else if beta<1 then result:=(x>=0) and (x<1)
  else result:=(x>=0) and (x<=1);
end;

function tBeta_distribution.density(x:extended):extended;
begin
  if (alpha=1) and (x=0) then result:=beta
  else if (beta=1) and (x=1) then result:=alpha
  else if (x>0) and (x<1) then
    result:=exp(const_part+(alpha-1)*ln(x)+(beta-1)*ln(1-x))
  else result:=0;
end;

function tBeta_distribution.cdf(x:extended):extended;
begin
  if in_support(x)
  then result:=incbeta(x,alpha,beta)
  else result:=inherited cdf(x);
end;

function tBeta_distribution.random_value:extended;
var u1,v,y:extended;
begin   {Devroye: Non-Uniform Random Variate Generation, p.438}
  repeat
    u1:=randomdraw; v:=ln(u1/(1-u1))/random_const1; y:=alpha*exp(v);
  until
    alpha_plus_beta*ln(alpha_plus_beta/(beta+y))+random_const2*v-ln(4)>=ln(sqr(u1)*randomdraw);
  result:=y/(beta+y);
end;

function tBeta_distribution.Mean:extended;
begin
  result:=alpha/(alpha_plus_beta);
end;

function tBeta_distribution.Variance:extended;
begin
  result:=mean*beta/alpha_plus_beta/(alpha_plus_beta+1);
end;

{-------------------------------------------------------------------}
{  tGamma_distribution                                              }
{-------------------------------------------------------------------}

constructor tGamma_distribution.create(newalpha, newlambda:extended);
begin
  Fparametercount:=2;
  with Fparameters[1] do
  begin
      bound.lower:=0; bound.upper:=maxvalue;
      integertype:=false;
  end;
  with Fparameters[2] do
  begin
    bound.lower:=0; bound.upper:=maxvalue;
    integertype:=false;
  end;
  Setparameters(newalpha,newlambda,0);
end;

function tGamma_distribution.Setparameters(par1,par2,par3:extended):boolean;
var z:extended;
begin
  result:=inherited Setparameters(par1,par2,par3);
  alpha:=par1; lambda:=par2;
  const_part := alpha*ln(lambda)-log_gamma(alpha);
  support.lower := 0; support.upper:=maxvalue;
  if alpha>1 then
    begin  {needed for simulation: random_const1=a-; random_const2:=a+}
      z:=sqrt(2*alpha-1);
      random_const1:=(1-z)*sqrt(exp((alpha-1)*(ln(alpha-z)-ln(alpha-1))+z-1));
      random_const2:=(1+z)*sqrt(exp((alpha-1)*(ln(alpha+z)-ln(alpha-1))-z-1));
    end
  else if alpha<1 then
    begin
      random_const1:=(exp(1)+alpha)/exp(1); random_const2:=1/alpha;
    end;
end;

function tGamma_distribution.in_support(x:extended):boolean;
begin
  if alpha<1
  then result:=inherited in_support(x) and (x>support.lower)
  else result:=inherited in_support(x);
end;

function tGamma_distribution.density(x:extended):extended;
begin
  if in_support(x) and (x>0)
  then result:=exp(const_part+(alpha-1)*ln(x)-lambda*x)
  else
    if approx(1,alpha) and (x=0)
    then result:=lambda
    else result:=0;
end;

function tGamma_distribution.cdf(x:extended):extended;
begin
  if in_support(x)
  then
    if alpha<=500
    then result:=incgamma(x*lambda, alpha)
    else
      result:=(1+incgamma(sqr(power(lambda*x/alpha,1/3)-1+1/(9*alpha))*4.5*alpha,0.5))/2
  else result:=inherited cdf(x);
end;

function tGamma_distribution.random_value:extended;
var u,v,w,x,x1,x2:extended; accept,qac,qrc:boolean;
begin
  if alpha=1  {exponential distribution}
  then x:=-ln(randomdraw)
  else
    if alpha>1  {unimodal}
    then  {Devroye: Non-Uniform Random Variate Generation, p.203}
      begin
        repeat
          u:=randomdraw; v:=random_const1+(random_const2-random_const1)*randomdraw;
          x:=v/u; x1:=x+alpha-1; x2:=sqr(x);
          if x>=0 then
            begin
               qac:=sqr(x1)*(-2*sqr(u)+8*u-6)<=-x2*(x+x1);
               qrc:=x1*(2*sqr(u)-2)>=-u*x2;
            end
          else
            begin
              qac:=x1*(-2*sqr(u)+8*u-6)<=-x2;
              qrc:=(alpha-1)*(2*sqr(u)-2)>=-u*x2;
            end;
        until qac or (not qrc and (x+alpha-1>0) and (2*ln(u)+x<=(alpha-1)*ln(1+x/(alpha-1))));
        x:=x+alpha-1;
      end
    else
      begin  {Devroye: Non-Uniform Random Variate Generation, p.425}
        repeat
          u:=randomdraw; v:=random_const1*u; w:=randomdraw;
          if v<=1
          then
            begin
              x:=power(v,random_const2);
              accept:=ln(w)<=-x;
            end
          else
            begin
              x:=-ln(random_const2*(random_const1-v));
              accept:=w<=power(x,alpha-1);
            end;
        until accept;
      end;
  result:=x/lambda;
end;

function tGamma_distribution.Mean:extended;
begin
  result:=alpha/lambda;
end;

function tGamma_distribution.Variance:extended;
begin
  result:=alpha/sqr(lambda);
end;

{-------------------------------------------------------------------}
{  tChi_square_distribution                                         }
{-------------------------------------------------------------------}

constructor tChi_square_distribution.create(newdf:extended);
begin
  Inherited create(newdf/2,1/2);
end;

function tChi_square_distribution.Setparameters(par1,par2,par3:extended):boolean;
begin
  result:=inherited Setparameters(par1,par2,par3);
  df:=2*par1;
end;

{-------------------------------------------------------------------}
{  tExponential_distribution                                        }
{-------------------------------------------------------------------}

constructor tExponential_distribution.create(newtheta:extended);
begin
  Inherited create(1,newtheta);
end;

function tExponential_distribution.Setparameters(par1,par2,par3:extended):boolean;
begin
  result:=inherited Setparameters(par1,par2,par3);
  Fparametercount:=1; Fparameters[1].value:=par2;
  theta:=par2;
end;

function tExponential_distribution.quantile(var q:extended):extended;
begin
  if q<=0 then result:=support.lower
  else if q>=1 then result:=support.upper
  else result:=-ln(1-q)/theta;
end;

{-------------------------------------------------------------------}
{  tF_distribution                                                  }
{-------------------------------------------------------------------}

constructor tF_distribution.create(newdf1,newdf2:extended);
begin
  Fparametercount:=2;
  with Fparameters[1] do
    begin
      bound.lower:=1; bound.upper:=maxvalue;
      integertype:=true;
    end;
  with Fparameters[2] do
    begin
      bound.lower:=1; bound.upper:=maxvalue;
      integertype:=true;
    end;
  cbeta:=tBeta_distribution.create(newdf1/2,newdf2/2);
  cp:=tChi_square_distribution.create(newdf1);
  Setparameters(newdf1,newdf2,0);
end;

destructor tF_distribution.free;
begin
  cp.free;
  cbeta.free;
  inherited free;
end;

function tF_distribution.Setparameters(par1,par2,par3:extended):boolean;
begin
  result:=inherited Setparameters(par1,par2,par3);
  df1:=par1; df2:=par2;
  const_part:=log_gamma((df1+df2)/2)-log_gamma(df1/2)-log_gamma(df2/2)+df1/2*(ln(df1)-ln(df2));
  support.lower := 0; support.upper:=maxvalue;
  cbeta.Setparameters(df1/2,df2/2,0);
  if df2>1000
  then cp.Setparameters(df1,0,0)
  else
    if df1>1000
    then cp.Setparameters(df2,0,0);
end;

function tF_distribution.in_support(x:extended):boolean;
begin
  if df1=1
  then result:=inherited in_support(x) and (x>support.lower)
  else result:=inherited in_support(x);
end;

function tF_distribution.density(x:extended):extended;
begin
  if in_support(x) and (x>0)
  then result:=cBeta.density(df1*x/(df1*x+df2))*df1*df2/sqr(df1*x+df2)
  else
    if (df1=2) and (x=0)
    then result:=exp(const_part)
    else result:=0;
end;

function tF_distribution.cdf(x:extended):extended;
begin
  if in_support(x) then
    begin
      if df2>1000 then result:=cp.cdf(df1*x)
      else if df1>1000 then result:=1-cp.cdf(df2/x)
      else result:=cBeta.cdf(x*df1/(df2+x*df1));
    end
  else result:=inherited cdf(x);
end;

function tF_distribution.random_value:extended;
var x:extended;
begin
  x:=cbeta.random_value;
  result:=df2*x/(df1*(1-x));
end;

function tF_distribution.expectation_exists:boolean;
begin
  result:=inherited expectation_exists and (df2>2);
end;

function tF_distribution.variance_exists:boolean;
begin
  result:=inherited expectation_exists and (df2>4);
end;

function tF_distribution.Mean:extended;
begin
  Check(not expectation_exists,err_expectation);
  result:=df2/ (df2 - 2);
end;

function tF_distribution.Variance:extended;
begin
  Check(not variance_exists,err_variance);
  result:=2*Sqr(Mean)*(df1+df2-2)/df1/(df2-4);
end;

{-------------------------------------------------------------------}
{  tT_distribution                                                  }
{-------------------------------------------------------------------}

constructor tT_distribution.create(newdf:extended);
begin
  Fparametercount:=1;
  With Fparameters[1] do
    begin
      bound.lower:=0; bound.upper:=maxvalue;
      integertype:=false;
    end;
  N:=tNormal_distribution.create(0,1);
  G:=tGamma_distribution.create(newdf/2,1);
  Setparameters(newdf,0,0);
end;

destructor tT_distribution.free;
begin
  N.free; G.free;
  inherited free;
end;

function tT_distribution.Setparameters(par1,par2,par3:extended):boolean;
begin
  result:=inherited Setparameters(par1,par2,par3);
  df:=par1;
  const_part:=log_gamma((df+1)/2)-log_gamma(df/2)-ln(df*pi)/2;
  if df>=1
  then random_const:=sqrt(2*df)*power((df-1),(df-1)/4)/power((df+1),(df+1)/4)
  else
    begin
      random_const:=sqrt(2*df);
      G.Setparameters(df/2,1,0);
    end;
end;

function tT_distribution.density(x:extended):extended;
begin
  if in_support(x)
  then result:=exp(const_part+(df+1)/2*ln(df/(df+sqr(x))))
  else result:=0;
end;

function tT_distribution.cdf(x:extended):extended;
begin
  if in_support(x)
  then
    begin
      if df>1000
      then
        if x<=0 then result:=(1-incgamma(sqr(x)/2,0.5))/2
        else result:=(1+incgamma(sqr(x)/2,0.5))/2
      else
        if x>=0 then result:=1 - 0.5*incbeta(df/(df+Sqr(x)), df/2, 1/2)
        else result:= 0.5*incbeta(df/(df+Sqr(x)), df/2, 1/2);
    end
  else result:=inherited cdf(x);
end;

function tT_distribution.random_value:extended;
var x,u1,u2:extended;
begin {Devroye: Non-Uniform Random Variate Generation, p.201}
  if (df=1) or (df=3) then
    begin
      repeat
        u1:=2*randomdraw-1; u2:=2*randomdraw-1;
      until sqr(u1)+sqr(u2)<=1;
      if df=1 then result:=u2/u1 else result:=sqrt(3)*u2/(1+u1);
    end
  else if df>1 then
    begin
      repeat
        u1:=randomdraw;
        u2:=random_const*(2*randomdraw-1);
        x:=u2/u1;
      until sqr(x)<=df*(power(u1,-4/(df+1))-1);
      result:=x;
    end
  else {if 0 < df < 1; Devroye, p.445}
    result:=random_const*N.random_value/sqrt(G.random_value);
end;

function tT_distribution.expectation_exists:boolean;
begin
  result:=inherited expectation_exists and (df>1);
end;

function tT_distribution.variance_exists:boolean;
begin
  result:=inherited expectation_exists and (df>2);
end;

function tT_distribution.Mean:extended;
begin
  Check(not expectation_exists,err_expectation);
  result:=0;
end;

function tT_distribution.Variance:extended;
begin
  Check(not variance_exists,err_variance);
  result:=df/(df-2);
end;

{-------------------------------------------------------------------}
{  tCauchy_distribution                                             }
{-------------------------------------------------------------------}

constructor tCauchy_distribution.create(newalpha, newbeta:extended);
begin
  Fparametercount:=2;
  With Fparameters[1] do
    begin
      bound.lower:=-maxvalue; bound.upper:=maxvalue;
      integertype:=false;
    end;
  With Fparameters[2] do
    begin
      bound.lower:=0; bound.upper:=maxvalue;
      integertype:=false;
    end;
  Setparameters(newalpha,newbeta,0);
end;

function tCauchy_distribution.Setparameters(par1,par2,par3:extended):boolean;
begin
  result:=inherited Setparameters(par1,par2,par3);
  alpha:=par1; beta:=par2;
end;

function tCauchy_distribution.density(x:extended):extended;
begin
  if in_support(x)
  then result:=beta/pi/(sqr(beta)+sqr(x-alpha))
  else result:=0;
end;

function tCauchy_distribution.cdf(x:extended):extended;
begin
  if in_support(x)
  then result:=0.5+arctan((x-alpha)/beta)/pi
  else result:=inherited cdf(x);
end;

function tCauchy_distribution.random_value:extended;
var v1,v2:extended;
begin
  repeat
    v1:=-1+2*randomdraw; v2:=-1+2*randomdraw;
  until sqr(v1)+sqr(v2)<=1;
  result:=v1/v2;
end;

function tCauchy_distribution.expectation_exists:boolean;
begin
  result:=false;
end;

function tCauchy_distribution.variance_exists:boolean;
begin
  result:=false;
end;

{-------------------------------------------------------------------}
{  double exponential's method implementations:                     }
{-------------------------------------------------------------------}

constructor tDouble_exponential_distribution.create(newalpha, newbeta:extended);
begin
  Fparametercount:=2;
  With Fparameters[1] do
    begin
      bound.lower:=-maxvalue; bound.upper:=maxvalue;
      integertype:=false;
    end;
  With Fparameters[2] do
    begin
      bound.lower:=0; bound.upper:=maxvalue;
      integertype:=false;
    end;
  Setparameters(newalpha,newbeta,0);
end;

function tDouble_exponential_distribution.Setparameters(par1,par2,par3:extended):boolean;
begin
  result:=inherited Setparameters(par1,par2,par3);
  alpha:=par1; beta:=par2;
end;

function tDouble_exponential_distribution.density(x:extended):extended;
begin
  if in_support(x)
  then result:=exp(-abs(x-alpha)/beta)/beta/2
  else result:=0;
end;

function tDouble_exponential_distribution.cdf(x:extended):extended;
begin
  if in_support(x)
  then
    if x<=alpha
    then result:=exp((x-alpha)/beta)/2
    else result:=1-exp((alpha-x)/beta)/2
  else result:=inherited cdf(x);
end;

function tDouble_exponential_distribution.random_value:extended;
begin
  if randomdraw<=0.5
  then result:=alpha + beta*ln(randomdraw)
  else result:=alpha - beta*ln(randomdraw);
end;

function tDouble_exponential_distribution.mean:extended;
begin
  result:=alpha;
end;

function tDouble_exponential_distribution.variance:extended;
begin
  result:=2*sqr(beta);
end;

{-------------------------------------------------------------------}
{  tGumbel_distribution                                             }
{-------------------------------------------------------------------}

constructor tGumbel_distribution.create(newalpha, newbeta:extended);
begin
  Fparametercount:=2;
  With Fparameters[1] do
    begin
      bound.lower:=-maxvalue; bound.upper:=maxvalue;
      integertype:=false;
    end;
  With Fparameters[2] do
    begin
      bound.lower:=0; bound.upper:=maxvalue;
      integertype:=false;
    end;
  Setparameters(newalpha,newbeta,0);
end;

function tGumbel_distribution.Setparameters(par1,par2,par3:extended):boolean;
begin
  result:=inherited Setparameters(par1,par2,par3);
  alpha:=par1; beta:=par2;
end;

function tGumbel_distribution.density(x:extended):extended;
begin
  if in_support(x)
  then result:=exp((alpha-x)/beta-exp((alpha-x)/beta))/beta
  else result:=0;
end;

function tGumbel_distribution.cdf(x:extended):extended;
begin
  if in_support(x)
  then result:=exp(-exp((alpha-x)/beta))
  else result:=inherited cdf(x);
end;

function tGumbel_distribution.random_value:extended;
begin
  result:=alpha-beta*ln(-ln(randomdraw));
end;


function tGumbel_distribution.mean:extended;
begin
  result:=alpha+0.577216*beta;
end;

function tGumbel_distribution.variance:extended;
begin
  result:=sqr(beta*pi)/6;
end;

{-------------------------------------------------------------------}
{  tLogistic_distribution                                           }
{-------------------------------------------------------------------}

constructor tLogistic_distribution.create(newalpha, newbeta:extended);
begin
  Fparametercount:=2;
  With Fparameters[1] do
    begin
      bound.lower:=-maxvalue; bound.upper:=maxvalue;
      integertype:=false;
    end;
  With Fparameters[2] do
    begin
      bound.lower:=0; bound.upper:=maxvalue;
      integertype:=false;
    end;
  Setparameters(newalpha,newbeta,0);
end;

function tLogistic_distribution.Setparameters(par1,par2,par3:extended):boolean;
begin
  result:=inherited Setparameters(par1,par2,par3);
  alpha:=par1; beta:=par2;
end;

function tLogistic_distribution.density(x:extended):extended;
var y:extended;
begin
  if in_support(x)
  then
    begin
      y:=exp((x-alpha)/beta);
      result:=1/(2+y+1/y)/beta;
    end
  else result:=0;
end;

function tLogistic_distribution.cdf(x:extended):extended;
begin
  if in_support(x)
  then result:=1/(1+exp((alpha-x)/beta))
  else result:=inherited cdf(x);
end;

function tLogistic_distribution.random_value:extended;
var u:extended;
begin
  u:=randomdraw;
  result:=alpha+beta*ln(u/(1-u));
end;

function tLogistic_distribution.mean:extended;
begin
  result:=alpha;
end;

function tLogistic_distribution.variance:extended;
begin
  result:=sqr(beta*pi)/3;
end;

{-------------------------------------------------------------------}
{  tLognormal_distribution                                          }
{-------------------------------------------------------------------}

constructor tLognormal_distribution.create(newmu, newsigma2:extended);
begin
  Fparametercount:=2;
  With Fparameters[1] do
    begin
      bound.lower:=-maxvalue; bound.upper:=maxvalue;
      integertype:=false;
    end;
  With Fparameters[2] do
    begin
      bound.lower:=0; bound.upper:=maxvalue;
      integertype:=false;
    end;
  norm:=tNormal_distribution.create(newmu,newsigma2);
  Setparameters(newmu,newsigma2,0);
end;

destructor tLognormal_distribution.free;
begin
  norm.free;
  inherited free;
end;

function tLognormal_distribution.Setparameters(par1,par2,par3:extended):boolean;
begin
  result:=inherited Setparameters(par1,par2,par3);
  mu:=par1; sigma2:=par2; sigma:=sqrt(sigma2);
  support.lower:=0;
  const_part:= -(ln(2*pi)+ln(sigma2))/2;
  norm.Setparameters(mu,sigma2,0);
end;

function tLognormal_distribution.in_support(x:extended):boolean;
begin
  result:=inherited in_support(x) and (x>0);
end;

function tLognormal_distribution.density(x:extended):extended;
begin
  if in_support(x)
  then result:=exp(const_part-ln(x)-sqr(ln(x)-mu)/sigma2/2)
  else result:=0;
end;

function tLognormal_distribution.cdf(x:extended):extended;
begin
  if in_support(x)
  then result:=norm.cdf(ln(x))
  else result:=inherited cdf(x);
end;

function tLognormal_distribution.random_value:extended;
begin
  result:=exp(norm.random_value);
end;

function tLognormal_distribution.mean:extended;
begin
  result:=exp(mu+sigma2/2);
end;

function tLognormal_distribution.variance:extended;
begin
  result:=exp(2*mu+2*sigma2)-exp(2*mu+sigma2);
end;

{-------------------------------------------------------------------}
{  tPareto_distribution                                             }
{-------------------------------------------------------------------}

constructor tPareto_distribution.create(newalpha, newtheta:extended);
begin
  Fparametercount:=2;
  With Fparameters[1] do
    begin
      bound.lower:=0; bound.upper:=maxvalue;
      integertype:=false;
    end;
  With Fparameters[2] do
    begin
      bound.lower:=0; bound.upper:=maxvalue;
      integertype:=false;
    end;
  Setparameters(newalpha,newtheta,0);
end;

function tPareto_distribution.Setparameters(par1,par2,par3:extended):boolean;
begin
  result:=inherited Setparameters(par1,par2,par3);
  alpha:=par1; theta:=par2;
  support.lower:=alpha;
end;

function tPareto_distribution.density(x:extended):extended;
begin
  if in_support(x)
  then result:=theta/alpha*power(alpha/x,theta+1)
  else result:=0;
end;

function tPareto_distribution.cdf(x:extended):extended;
begin
  if in_support(x)
  then result:=1-power(alpha/x,theta)
  else result:=inherited cdf(x);
end;

function tPareto_distribution.random_value:extended;
begin
  result:=alpha*power(randomdraw,-1/theta);
end;

function tPareto_distribution.expectation_exists:boolean;
begin
  result:=inherited expectation_exists and (theta>1);
end;

function tPareto_distribution.variance_exists:boolean;
begin
  result:=inherited variance_exists and (theta>2);
end;

function tPareto_distribution.mean:extended;
begin
  Check(not expectation_exists,err_expectation);
  result:=alpha*theta/(theta-1);
end;

function tPareto_distribution.variance:extended;
begin
  Check(not variance_exists,err_variance);
  result:=mean*alpha/((theta-1)*(theta-2));
end;

{-------------------------------------------------------------------}
{  tUniform_distribution                                            }
{-------------------------------------------------------------------}

constructor tUniform_distribution.create(newalpha, newbeta:extended);
begin
  Fparametercount:=2;
  With Fparameters[1] do
    begin
      bound.lower:=-maxvalue; bound.upper:=newbeta;
      integertype:=false;
    end;
  With Fparameters[2] do
    begin
      bound.lower:=newalpha; bound.upper:=maxvalue;
      integertype:=false;
    end;
  Setparameters(newalpha,newbeta,0)
end;

function tUniform_distribution.Setparameters(par1,par2,par3:extended):boolean;
begin
  Fparameters[1].bound.upper:=par2;
  Fparameters[2].bound.lower:=par1;
  result:=inherited Setparameters(par1,par2,par3);
  alpha:=par1; beta:=par2;
  support.lower:=alpha; support.upper:=beta;
end;

function tUniform_distribution.density(x:extended):extended;
begin
  if in_support(x)
  then result:=1/(beta-alpha)
  else result:=0;
end;

function tUniform_distribution.cdf(x:extended):extended;
begin
  if in_support(x)
  then result:=(x-alpha)/(beta-alpha)
  else result:=inherited cdf(x);
end;

function tUniform_distribution.random_value:extended;
begin
  result:=alpha+(beta-alpha)*randomdraw;
end;

function tUniform_distribution.mean:extended;
begin
  result:=(alpha+beta)/2;
end;

function tUniform_distribution.variance:extended;
begin
  result:=sqr(beta-alpha)/12;
end;

{-------------------------------------------------------------------}
{  tWeibull_distribution                                            }
{-------------------------------------------------------------------}

constructor tWeibull_distribution.create(newa, newb:extended);
begin
  Fparametercount:=2;
  With Fparameters[1] do
    begin
      bound.lower:=0; bound.upper:=maxvalue;
      integertype:=false;
    end;
  With Fparameters[2] do
    begin
      bound.lower:=0; bound.upper:=maxvalue;
      integertype:=false;
    end;
  Setparameters(newa,newb,0);
end;

function tWeibull_distribution.Setparameters(par1,par2,par3:extended):boolean;
begin
  result:=inherited Setparameters(par1,par2,par3);
  a:=par1; b:=par2; b1:=1/b;
  support.lower:=0;
end;

function tWeibull_distribution.in_support(x:extended):boolean;
begin
  result:=inherited in_support(x);
  if b<=1
  then result:=result and (x>0);
end;

function tWeibull_distribution.density(x:extended):extended;
begin
  if in_support(x) and (x>0)
  then result:=a*b*exp((b-1)*ln(x)-a*power(x,b))
  else result:=0;
end;

function tWeibull_distribution.cdf(x:extended):extended;
begin
  if in_support(x)
  then result:=1-exp(-a*power(x,b))
  else result:=inherited cdf(x);
end;

function tWeibull_distribution.random_value:extended;
begin
  result:=power(-ln(randomdraw),b1);
end;

function tWeibull_distribution.mean:extended;
begin
  result:=exp(-b1*ln(a)+log_gamma(1+b1));
end;

function tWeibull_distribution.variance:extended;
begin
  result:=power(a,-2*b1)*(exp(log_gamma(1+2*b1))-sqr(exp(log_gamma(1+b1))));
end;

{-------------------------------------------------------------------}
{  tNon_Central_Chi_square_distribution                             }
{-------------------------------------------------------------------}

constructor tNon_Central_Chi_square_distribution.create(newdf,newlambda:extended);
begin
  Fparametercount:=2;
  With Fparameters[1] do
    begin
      bound.lower:=1; bound.upper:=maxvalue;
      integertype:=true;
    end;
  With Fparameters[2] do
    begin
      bound.lower:=0; bound.upper:=maxvalue;
      integertype:=false;
    end;
  stand_norm:=tNormal_distribution.create(0,1);
  central_chi2:=tChi_square_distribution.create(newdf);
  c_chi2:=tchi_square_distribution.create(newdf);
  Setparameters(newdf,newlambda,0);
  support.lower := 0;
end;

destructor tNon_Central_Chi_square_distribution.free;
begin
  stand_norm.free;
  central_chi2.free;
  c_chi2.free;
  inherited free;
end;

function tNon_Central_Chi_square_distribution.Setparameters(par1,par2,par3:extended):boolean;
begin
  result:=inherited Setparameters(par1,par2,par3);
  df:=par1; lambda:=par2;
  central_chi2.Setparameters(df,0,0);
end;

function tNon_Central_Chi_square_distribution.density(x:extended):extended;
var j,k:integer; c,pk,p,fk,f,contribution:extended;
begin
  if x<=0 then result:=0
  else
    if approx(lambda,0)
    then result:=central_chi2.density(x)
  else
  begin
  c:=lambda/2;
  k:=round(c);
  pk:=exp(-c+k*ln(c)-log_gamma(k+1));
  fk:=exp(-log_gamma(df/2+k)-(df/2+k)*ln(2)-x/2+(df/2+k-1)*ln(x));
  j:=k; p:=pk; f:=fk; result:=pk*fk; contribution:=result;
  while (j>0) and (contribution>eps) do
    begin
      dec(j);
      p:=p*(j+1)/c;
      f:=f*(df+2*j)/x;
      contribution:=p*f;
      result:=result+contribution;
    end;
  j:=k; p:=pk; f:=fk;
  repeat
    inc(j);
    p:=p*c/j;
    f:=f*x/(df+2*j-2);
    contribution:=p*f;
    result:=result+contribution;
  until contribution<eps;
end;
end;

function tNon_Central_Chi_square_distribution.cdf(x:extended):extended;
var c,pk,Fk,Gk,p,F,G,sump:extended; j,k,m:integer;
begin
if (x>0)
then
  if approx(lambda,0) then result:=central_chi2.cdf(x)
  else
    begin
      c:=lambda/2;
      m:=floor(c-8*sqrt(c)); if m<0 then m:=0;
      k:=round(c); if k<m then k:=m;
      c_chi2.Setparameters(df+2*k,0,0);
      Fk:=c_chi2.cdf(x);
      pk:=exp(-c+k*ln(c)-log_gamma(k+1));
      Gk:=exp(-x/2+(df/2+k)*ln(x/2)-log_gamma(df/2+k+1));
      result:=pk*Fk;
      p:=pk; F:=Fk; G:=Gk; sump:=1-p;
      for j:=k-1 downto m do
        begin
        p:=p*(j+1)/c; sump:=sump-p;
        G:=G*(df/2+j+1)/x*2;
        F:=F+G;
        result:=result+p*F;
        end;
      p:=pk; F:=Fk; G:=Gk; j:=k;
      repeat
        inc(j);
        p:=p*c/j; sump:=sump-p;
        F:=F-G;
        G:=G*x/2/(df/2+j);
        result:=result+p*F;
      until (sump*(F-G)<eps);
    end
  else result:=inherited cdf(x);
end;

function tNon_Central_Chi_square_distribution.mean:extended;
begin
    result:=df+lambda;
end;

function tNon_Central_Chi_square_distribution.variance:extended;
begin
    result:=2*(df+2*lambda);
end;


{-------------------------------------------------------------------}
{  tNon_Central_Beta_distribution                                   }
{-------------------------------------------------------------------}

{ For the methods `density' and `cdf' Rich Strauss <y8res@ttacs1.ttu.edu>
  pointed me to the following papers:

Norton, V.  1983.  A simple algorithm for computing the con-central F
distribution.  Appl. Stat. 32:84-85.

Lenth, R.V.  1987.  Algorithm AS 226.  Computing non-central beta
probabilities.  Appl. Stat. 36:241-244.

Frick, H.  1990.  Algorithm AS R84.  A remark on Algorithm AS 226:
computing non-central beta probabilities.  Appl. Stat. 39:311-312. }

constructor tNon_Central_Beta_distribution.create(newalpha,newbeta,newlambda:extended);
begin
  Fparametercount:=3;
  With Fparameters[1] do
    begin
      bound.lower:=0; bound.upper:=maxvalue;
      integertype:=false;
    end;
  With Fparameters[2] do
    begin
      bound.lower:=0; bound.upper:=maxvalue;
      integertype:=false;
    end;
  With Fparameters[3] do
    begin
      bound.lower:=0; bound.upper:=maxvalue;
      integertype:=false;
    end;
  Setparameters(newalpha,newbeta,newlambda);
  support.lower := 0; support.upper :=1;
end;

function tNon_central_beta_distribution.Setparameters(par1,par2,par3:extended):boolean;
begin
  result:=inherited Setparameters(par1,par2,par3);
  alpha:=par1; beta:=par2; lambda:=par3;
end;

function tNon_Central_Beta_distribution.in_support(x:extended):boolean;
begin
  if (alpha<1) and (beta<1) then result:=(x>0) and (x<1)
  else if alpha<1 then result:=(x>0) and (x<=1)
  else if beta<1 then result:=(x>=0) and (x<1)
  else result:=(x>=0) and (x<=1);
end;

function tNon_Central_Beta_distribution.density(x:extended):extended;
var m,k,j:integer; c,qk,sumq,fk,qj,fj:extended;
begin
  if in_support(x)
  then
    begin
      c := lambda/2;
      m:=floor(c-8*sqrt(c)); if m<0 then m:=0;
      k:=round(c); if k<m then k:=m;
      if k=0
      then qk:=exp(-c)
      else qk:=exp(-c+k*ln(c)-log_gamma(k+1));
      sumq:=1-qk;
      fk:=exp((alpha+k-1)*Ln(x)+(beta-1)*Ln(1-x)
          +log_gamma(alpha+k+beta)-log_gamma(alpha+k)-log_gamma(beta));
      result:=qk*fk;   {contribution for j=k}
      j:=k;
      qj:=qk; fj:=fk;
      while j>m do {contributions for m <= j < k}
        begin
          dec(j);
          fj:=fj*(alpha+j)/(alpha+beta+j)/x;
          qj:=qj*(j+1)/c; sumq:=sumq-qj;
          result:=result+fj*qj;
        end;
      j:=k;
      qj:=qk; fj:=fk;
      repeat       {contributions for j > k}
        inc(j);
        fj:=fj*x*(alpha+beta+j-1)/(alpha+j-1);
        qj:=qj*c/j; sumq:=sumq-qj;
        result:=result+fj*qj;
      until (fj*sumq<eps);
    end
 else result:=0;
end;

function tNon_Central_Beta_distribution.cdf(x:extended):extended;
var j,k,m:integer; c,qk,sumq,Gk,Ik,qj,Gj,Ij:extended;
begin
  if (x>0) and (x<1)
  then
    begin
      c := lambda/2;
      m:=floor(c-8*sqrt(c)); if m<0 then m:=0;
      k:=round(c); if k<m then k:=m;
      if k=0
      then qk:=exp(-c)
      else qk:=exp(-c+k*ln(c)-log_gamma(k+1));
      sumq:=1-qk;
      Gk:=exp((alpha+k)*Ln(x)+beta*Ln(1-x)
          +log_gamma(alpha+k+beta)-log_gamma(alpha+k+1)-log_gamma(beta));
      Ik:=incbeta(x,alpha+k,beta);
      result:=qk*Ik;   {contribution for j=k}
      qj:=qk; Gj:=Gk; Ij:=Ik;
      for j:=k-1 downto m do   {contributions for j < k}
        begin
          Gj:=Gj*(alpha+j+1)/(alpha+beta+j)/x;
          Ij:=Ij+Gj;
          qj:=qj*(j+1)/c; sumq:=sumq-qj;
          result:=result+Ij*qj;
        end;
      j:=k; qj:=qk; Gj:=Gk; Ij:=Ik;
      repeat       {contributions for j > k}
        j:=j+1;
        Ij:=Ij-Gj;
        Gj:=Gj*x*(alpha+beta+j-1)/(alpha+j);
        qj:=qj*c/j; sumq:=sumq-qj;
        result:=result+Ij*qj;
      until ((Ij-Gj)*sumq<eps);
    end
  else result:=inherited cdf(x);
end;

function tNon_Central_Beta_distribution.mean:extended;
begin
  result:=alpha/(alpha+beta);
end;

function tNon_Central_Beta_distribution.variance:extended;
begin
  result:= mean*beta/(alpha+beta)/(alpha+beta+1);
end;

{-------------------------------------------------------------------}
{  tNon-central_F_distribution                                      }
{-------------------------------------------------------------------}

constructor tNon_central_F_distribution.create(newdf1,newdf2,newlambda:extended);
begin
  Fparametercount:=3;
  With Fparameters[1] do
    begin
      bound.lower:=0; bound.upper:=maxvalue;
      integertype:=false;
    end;
  With Fparameters[2] do
    begin
      bound.lower:=0; bound.upper:=maxvalue;
      integertype:=false;
    end;
  With Fparameters[3] do
    begin
      bound.lower:=0; bound.upper:=maxvalue;
      integertype:=false;
    end;
  ncBeta:=tNon_central_Beta_distribution.Create(newdf1/2,newdf2/2,newlambda);
  Setparameters(newdf1,newdf2,newlambda);
  support.lower := 0;
end;

destructor tNon_central_F_distribution.free;
begin
  ncBeta.free;
  inherited free;
end;

function tNon_central_F_distribution.Setparameters(par1,par2,par3:extended):boolean;
begin
  result:=inherited Setparameters(par1,par2,par3);
  df1:=par1; df2:=par2; lambda:=par3;
  ncBeta.Setparameters(df1/2,df2/2,lambda);
end;

function tNon_central_F_distribution.in_support(x:extended):boolean;
begin
  result:=inherited in_support(x);
  if df1=1 then result:=result and (x>support.lower);
end;

function tNon_central_F_distribution.density(x:extended):extended;
begin
  if in_support(x) then
  result:=ncBeta.density(df1*x/(df1*x+df2))*df1*df2/sqr(df1*x+df2)
  else result:=0;
end;

function tNon_central_F_distribution.cdf(x:extended):extended;
begin
  if in_support(x)
  then result:=ncBeta.cdf(x*df1/(df2+x*df1))
  else result:=inherited cdf(x);
end;

function tNon_central_F_distribution.expectation_exists:boolean;
begin
  result := df2>2;
end;

function tNon_central_F_distribution.variance_exists:boolean;
begin
  result := df2>4;
end;

function tNon_central_F_distribution.mean:extended;
begin
  result:= df2*(df1+lambda)/(df1*(df2 - 2));
end;

function tNon_central_F_distribution.variance:extended;
begin
  result:=2*Sqr(df2/df1)*(sqr(df1+lambda)+(df1+2*lambda)*(df2-2))/(Sqr(df2-2)*(df2-4));
end;

{-------------------------------------------------------------------}
{  tNon_central_t_distribution                                      }
{-------------------------------------------------------------------}

constructor tNon_central_t_distribution.create(newdf,newdelta:extended);
begin
  Fparametercount:=2;
  With Fparameters[1] do
    begin
      bound.lower:=0; bound.upper:=maxvalue;
      integertype:=false;
    end;
  With Fparameters[2] do
    begin
      bound.lower:=-maxvalue; bound.upper:=maxvalue;
      integertype:=false;
    end;
  central_t:=tT_distribution.create(newdf);
  Setparameters(newdf,newdelta,0);
end;

destructor tNon_central_t_distribution.free;
begin
  central_t.free;
  inherited free;
end;

function tNon_central_t_distribution.Setparameters(par1,par2,par3:extended):boolean;
begin
  result:=inherited Setparameters(par1,par2,par3);
  df:=par1; delta:=par2;
  central_t.Setparameters(df,0,0);
end;

function tNon_central_t_distribution.density(x:extended):extended;
var j,k,m,negx:integer; negdelta:boolean;
    del,a,b,c,pk,qk,f1k,f2k,p,q,f1,f2:extended;
begin
  if approx(0,delta) then result:=central_t.density(x)
  else
    if in_support(x)
    then
    begin
    negdelta:=delta<0;
    if negdelta
    then
      begin
        x:=-x;
        del:=-delta;
      end
    else del:=delta;
    a:=sqr(delta)/2;
    m:=floor(a-8*sqrt(a)); if m<0 then m:=0;
    k:=round(a); if k<m then k:=m;
    if x=0 then k:=0;
    b:=df/2;
    c:=sqr(x)/(sqr(x)+df);
    pk:=exp(-a+k*ln(a)-log_gamma(k+1));
    if x<0 then negx:=-1 else negx:=1;
    qk:=del*exp(-a+k*ln(a)-0.5*ln(2)-log_gamma(k+1.5));
    if x=0
    then
      if k>0 then f1k:=0
      else f1k:=exp(log_gamma(0.5+b)-log_gamma(0.5)-log_gamma(b)-0.5*ln(df))
    else
      f1k:=exp(log_gamma(k+0.5+b)-log_gamma(k+0.5)-log_gamma(b))
         *power(c,k-0.5)*power(1-c,b-1)*df/sqr(sqr(x)+df)*x;
    if x=0
    then f2k:=0
    else
      f2k:=exp(log_gamma(k+1+b)-log_gamma(k+1)-log_gamma(b))
         *power(c,k)*power(1-c,b-1)*df/sqr(sqr(x)+df)*x;
    p:=pk; q:=qk; f1:=f1k; f2:=f2k;
    result:=negx*p*f1+q*f2;
    for j:=k-1 downto m do
      begin
        p:=p*(j+1)/a; q:=q*(j+1.5)/a;
        if x<>0
          then
            begin
        f1:=f1*(sqr(x)+df)/sqr(x)*(j+0.5)/(j+0.5+b);
        f2:=f2*(sqr(x)+df)/sqr(x)*(j+1)/(j+1+b);
            end;
        result:=result+negx*p*f1+q*f2;
      end;
    p:=pk; q:=qk; f1:=f1k; f2:=f2k; j:=k;
    repeat
      inc(j);
      p:=p*a/j; q:=q*a/(j+0.5);
      f1:=f1*sqr(x)/(sqr(x)+df)*(j-0.5+b)/(j-0.5);
      f2:=f2*sqr(x)/(sqr(x)+df)*(j+b)/j;
      result:=result+negx*p*f1+q*f2;
    until (abs(q*f2)<eps);
  end
  else result:=0;
end;

function tNon_central_t_distribution.cdf(x:extended):extended;
{ Based on Algorithm AS 243:
  Russell V. Lenth:
  Cumulative Distribution Function of the Non-central t Distribution,
  Applied Statistics (1989), Vol. 38, No.1, pp. 185 -189. }
var a,b,p,q,sump,I1,I2,G1,G2,del,pk,qk,I1k,I2k,G1k,G2k:extended;
    j,m,k,negx:integer; negdelta:boolean;
begin
  if approx(delta,0) then result:=central_t.cdf(x)
  else if in_support(x)
  then
  begin
    negdelta:=delta<0;
    if negdelta
    then
      begin
        x:=-x;
        del:=-delta;
      end
    else del:=delta;
    a:=sqr(del)/2;
    m:=floor(a-8*sqrt(a)); if m<0 then m:=0;
    k:=round(a); if k<m then k:=m;
    b:=df/2;
    pk:=0.5*exp(-a+k*ln(a)-log_gamma(k+1));
    if x<0 then negx:=-1 else negx:=1;
    qk:=0.5*del*exp(-a+k*ln(a)-0.5*ln(2)-log_gamma(k+1.5));
    x:=sqr(x)/(sqr(x)+df);
    I1k:=incbeta(x,k+0.5,b); I2k:=incbeta(x,k+1,b);
    G1k:=exp(log_gamma(k+0.5+b)-log_gamma(k+1.5)-log_gamma(b))*power(x,k+0.5)*power(1-x,b);
    G2k:=exp(log_gamma(k+1+b)-log_gamma(k+2)-log_gamma(b))*power(x,k+1)*power(1-x,b);
    p:=pk; q:=qk; I1:=I1k; I2:=I2k; G1:=G1k; G2:=G2k; sump:=0.5-p;
    result:=0.5*(1-incgamma(sqr(del)/2,0.5))+negx*p*I1+q*I2;
    for j:=k-1 downto m do
      begin
        p:=p*(j+1)/a; q:=q*(j+1.5)/a; sump:=sump-p;
        if x>0 then G1:=G1*(j+1.5)/(j+0.5+b)/x;
        if x>0 then G2:=G2*(j+2)/(j+1+b)/x;
        I1:=I1+G1; I2:=I2+G2;
        result:=result+negx*p*I1+q*I2;
      end;
    p:=pk; q:=qk; I1:=I1k; I2:=I2k; G1:=G1k; G2:=G2k; j:=k;
    repeat
      inc(j);
      p:=p*a/j; q:=q*a/(j+0.5); sump:=sump-p;
      I1:=I1-G1; I2:=I2-G2;
      result:=result+negx*p*I1+q*I2;
      G1:=G1*x*(j-0.5+b)/(j+0.5);
      G2:=G2*x*(j+b)/(j+1);
    until (sump*(I1-G1)<0.5*eps);
    if negdelta then result:=1-result;
  end
  else result:=inherited cdf(x);
end;

function tNon_central_t_distribution.expectation_exists:boolean;
begin
    result:=df>1;
end;

function tNon_central_t_distribution.variance_exists:boolean;
begin
    variance_exists:= df>2;
end;

function tNon_central_t_distribution.mean:extended;
begin
    result:=delta*sqrt(df/2)*exp(log_gamma((df-1)/2)-log_gamma(df/2));
end;

function tNon_central_t_distribution.variance:extended;
begin
    result:=df/(df-2)*(1+sqr(delta))-sqr(mean);
end;

{-------------------------------------------------------------------}
{  tDiscrete_distribution's method implementations:                  }
{-------------------------------------------------------------------}

procedure tDiscrete_distribution.fill_constants;
begin
  mode:=round(mean);
  while probability(mode-1)>probability(mode) do mode:=mode-1;
  while probability(mode+1)>probability(mode) do mode:=mode+1;
  maxprob:=probability(mode);
  random_const:=power(3*(variance+sqr(mode-mean))*sqr(maxprob),1/3);
end;

function tDiscrete_distribution.in_support(x:extended):boolean;
begin
  result:=inherited in_support(x) and nearint(x);
end;

function tDiscrete_distribution.prob_between(a,b:extended):extended; {P(a<=X<=b)}
var i,ia,ib:integer; prob,grootste,som:extended;
begin
  if a<support.lower then a:=support.lower;
  if in_support(a) then ia:=round(a) else ia:=floor(a+1);
  if b>support.upper then b:=support.upper;
  if in_support(b) then ib:=round(b) else ib:=floor(b);
  if ia>ib then prob_between:=0
  else
    if ia-support.lower+support.upper-ib<ib-ia+1 then
      begin
        if b>0.1*maxlongint then
        result:=1-prob_between(support.lower,ia-1) else
        result:=1-prob_between(support.lower,ia-1)-prob_between(ib+1,support.upper)
      end
    else
      begin
        i:=ia; prob:=probability(i); som:=prob; grootste:=prob;
        while (i<ib) and ((prob>eps) or (prob>=grootste)) do
          begin
            inc(i);
            if factor(i)>0 then prob:=prob * factor(i)
            else prob:=probability(i);
            som:=som+prob;
            if prob>grootste then grootste:=prob;
          end;
        if som<0 then result:=0 else
        if som>1 then result:=1 else
        result:=som;
      end;
end;

function tDiscrete_distribution.factor(i:integer):extended;
begin
  result:=0;
end;

function tDiscrete_distribution.quantile(var q:extended):extended;
Var x:integer; prob, som, fact: extended;
begin
  if (not Checkval(q,err_val_l,0) or approx(0,q)) then quantile:=support.lower
  else if (not Checkval(q,err_val_g,1) or approx(1,q)) then quantile:=support.upper
  else
    begin
      if expectation_exists then x:=round(mean)
      else x:=round(support.lower);
      if variance_exists then
        if q>0.5 then inc(x,round(sqrt(-variance*ln(2*(1-q)))))
        else dec(x,round(sqrt(-variance*ln(2*q))))
      else inc(x);
      if x>support.upper then x:=round(support.upper);
      if x<support.lower then x:=round(support.lower);
      som:=cdf(x); prob:=probability(x);
      while som>q+eps do
        begin
          som := som - prob; fact:=factor(x); dec(x);
          if fact > 0 then prob:=prob/fact else prob:=probability(x);
        end;
      while som<q-eps do
        begin
          inc(x); fact:=factor(x);
          if fact=0 then prob:=probability(x) else prob:=prob*fact;
          som := som + prob;
        end;
      q := som; result:= x;
    end;
end;

function tDiscrete_distribution.random_value:extended;
var u,v,w,x,t:extended;
begin
  repeat
    u:=randomdraw; v:=2*randomdraw-1; w:=randomdraw;
    if u<random_const/(3*random_const+maxprob)
    then
      begin
        if v>0
        then x:=round(mode+(0.5+random_const/(maxprob*sqrt(abs(v)))))
        else x:=round(mode-(0.5+random_const/(maxprob*sqrt(abs(v)))));
        t:=w*maxprob*abs(v)*sqrt(abs(v));
      end
    else
      begin
        x:=round(mode+(0.5+random_const/maxprob)*v);
        t:=w*maxprob;
      end;
  until t<=probability(x);
  result:=x;
end;

{-------------------------------------------------------------------}
{  tBinomial_distribution                                            }
{-------------------------------------------------------------------}

constructor tBinomial_distribution.create(newn, newp:extended);
begin
  Fparametercount:=2;
  With Fparameters[1] do
    begin
      bound.lower:=1; bound.upper:=maxlongint;
      integertype:=true;
    end;
  With Fparameters[2] do
    begin
      bound.lower:=0; bound.upper:=1;
      integertype:=false;
    end;
  Setparameters(newn,newp,0);
end;

function tBinomial_distribution.Setparameters(par1,par2,par3:extended):boolean;
begin
  result:=inherited Setparameters(par1,par2,par3);
  n:=round(par1); p:=par2;
  support.lower:=0; support.upper:=n;
  fill_constants;
end;

function tBinomial_distribution.factor(i:integer):extended;
begin
  if i>0 then result:=(n-i+1)*p/i/(1-p) else result:=0;
end;

function tBinomial_distribution.probability(x:extended): extended;
var ix:integer;
begin
  result:=0;
  if in_support(x) then
  begin
    ix:=round(x);
    if nearint(p) then
      case round(p) of
        0: if ix=0 then result:=1;
        1: if ix=n then result:=1;
      end
    else result:=exp(log_choose(n,ix)+ix*ln(p)+(n-ix)*ln(1-p));
  end;
end;

function tBinomial_distribution.cdf(x:extended):extended;
begin
  if (x>=support.lower) and (x<support.upper)
  then
    begin
      if nearint(x) then x:=round(x) else x:=floor(x);
      result:=1-incbeta(p,x+1,n-x);
    end
  else result:=inherited cdf(x);
end;

function tBinomial_distribution.mean:extended;
begin
  result:=n*p;
end;

function tBinomial_distribution.variance:extended;
begin
  result:=mean*(1-p);
end;

{-------------------------------------------------------------------}
{  tBernoulli_distribution                                          }
{-------------------------------------------------------------------}

constructor tBernoulli_distribution.create(newp:extended);
begin
  Inherited create(1,p);
end;

function tBernoulli_distribution.quantile(var q:extended):extended;
begin
  if (not Checkval(q,err_val_l,0) or approx(0,q)) then quantile:=support.lower
  else if (not Checkval(q,err_val_g,1) or approx(1,q)) then quantile:=support.upper
  else
    if q<=1-p then
      begin
        q:=1-p; result:=0;
      end
    else
      begin
        q:=1; result:=1;
      end;
end;

function tBernoulli_distribution.random_value:extended;
begin
  if randomdraw<p then result:=1 else result:=0;
end;

{-------------------------------------------------------------------}
{  tHypergeometric_distribution                                     }
{-------------------------------------------------------------------}

constructor tHypergeometric_distribution.create(newpopul_N, newN1, newsample_n:extended);
begin
  Fparametercount:=3;
  With Fparameters[1] do
    begin
      bound.lower:=2; bound.upper:=maxvalue;
      integertype:=true;
    end;
  With Fparameters[2] do
    begin
      bound.lower:=1; bound.upper:=newpopul_N;
      integertype:=true;
     end;
  With Fparameters[3] do
    begin
      bound.lower:=1; bound.upper:=newpopul_N;
      integertype:=true;
    end;
  Setparameters(newpopul_N, newN1, newsample_n);
end;

function tHypergeometric_distribution.Setparameters(par1,par2,par3:extended):boolean;
begin
  Fparameters[2].bound.upper:=par1;
  Fparameters[3].bound.lower:=par1;
  result:=inherited Setparameters(par1,par2,par3);
  popul_N:=round(par1); N1:=round(par2); N2:=popul_N-N1;
  if sample_n<=N2 then support.lower:=0 else support.lower:=sample_n-N2;
  if sample_n<=N1 then support.upper:=sample_n else support.upper:=N1;
  fill_constants;
end;

function tHypergeometric_distribution.probability(x:extended): extended;
var i,ix:integer;
begin
  if in_support(x) then
  begin
    ix:=round(x);
    result:=round(exp(log_choose(sample_n,ix)));
    for i:=0 to ix-1 do result:=result*(N1-i)/(popul_N-i);
    for i:=0 to sample_n-ix-1 do result:=result*(N2-i)/(popul_N-ix-i);
  end
  else result:=0;
end;

function tHypergeometric_distribution.factor(i:integer):extended;
begin
  if (i=0) or (N2+i=sample_n) then result:=0
  else result:=(N1-i+1)/i*(sample_n-i+1)/(N2-sample_n+i);
end;

function tHypergeometric_distribution.mean:extended;
begin
  result:=sample_n*N1/popul_N;
end;

function tHypergeometric_distribution.variance:extended;
begin
  result:=mean*N2*(popul_N-sample_n)/popul_N/(popul_N-1);
end;

{-------------------------------------------------------------------}
{  tPoisson_distribution                                            }
{-------------------------------------------------------------------}

constructor tPoisson_distribution.create(newlambda:extended);
begin
  Fparametercount:=1;
  With Fparameters[1] do
    begin
      bound.lower:=0; bound.upper:=maxvalue;
      integertype:=false;
    end;
  Setparameters(newlambda,0,0);
  support.lower := 0;
end;

function tPoisson_distribution.Setparameters(par1,par2,par3:extended):boolean;
begin
  result:=inherited Setparameters(par1,par2,par3);
  lambda:=par1;
  fill_constants;
end;

function tPoisson_distribution.probability(x:extended):extended;
var ix:integer;
begin
  if in_support(x)
  then
    begin
      ix:=round(x);
      if approx(0,lambda) then
        begin
          if ix=0 then result:=1 else result:=0;
        end
      else result:=exp(ix*ln(lambda)-lambda-log_gamma(ix+1));
    end
  else result:=0;
end;

function tPoisson_distribution.cdf(x:extended):extended;
begin
  if nearint(x) then x:=round(x) else x:=floor(x);
  if x<0
  then result:=0
  else result:=1-incgamma(lambda,x+1);
end;

function tPoisson_distribution.mean:extended;
begin
  result:=lambda;
end;

function tPoisson_distribution.variance:extended;
begin
  result:=lambda;
end;

{-------------------------------------------------------------------}
{  tDiscrete_uniform_distribution                                   }
{-------------------------------------------------------------------}

constructor tDiscrete_uniform_distribution.create(newN:extended);
begin
  Fparametercount:=1;
  With Fparameters[1] do
    begin
      bound.lower:=2; bound.upper:=maxvalue;
      integertype:=true;
    end;
  Setparameters(newN,0,0);
  support.lower:=1; support.upper:=N;
end;

function tDiscrete_uniform_distribution.Setparameters(par1,par2,par3:extended):boolean;
begin
  result:=inherited Setparameters(par1,par2,par3);
  N:=round(par1);
  fill_constants;
end;

function tDiscrete_uniform_distribution.probability(x:extended):extended;
begin
  if in_support(x) then result:=1/N else result:=0;
end;

function tDiscrete_uniform_distribution.random_value:extended;
begin
  result:=round(N*randomdraw+0.5);
end;

function tDiscrete_uniform_distribution.mean:extended;
begin
  result:=(N+1)/2;
end;

function tDiscrete_uniform_distribution.variance:extended;
begin
  result:=(sqr(N)-1)/12;
end;

{-------------------------------------------------------------------}
{  tNegative_binomial_distribution                                  }
{-------------------------------------------------------------------}

constructor tNegative_binomial_distribution.create(newr,newp:extended);
begin
  Fparametercount:=2;
  With Fparameters[1] do
     begin
      bound.lower:=0; bound.upper:=maxlongint;
      integertype:=true;
    end;
  with Fparameters[2] do
    begin
      bound.lower:=0; bound.upper:=1;
      integertype:=false;
    end;
  Setparameters(newr,newp,0);
  support.lower:=0;
end;

function tNegative_binomial_distribution.Setparameters(par1,par2,par3:extended):boolean;
begin
  result:=inherited Setparameters(par1,par2,par3);
  r:=round(par1); p:=par2;
  fill_constants;
end;

function tNegative_binomial_distribution.factor(i:integer):extended;
begin
  if i>0 then result:=(r+i-1)*(1-p)/i else result:=0;
end;

function tNegative_binomial_distribution.probability(x:extended): extended;
var ix:integer;
begin
  if in_support(x) then
  begin
    ix:=round(x);
    if approx(0,p)
    then result:=0
    else
      if ix=0
      then result:=exp(r*ln(p))
      else
        if approx(1,p)
        then result:=0
        else result:=exp(log_choose(r+ix-1,ix)+r*ln(p)+ix*ln(1-p));
  end
  else result:=0;
end;

function tNegative_binomial_distribution.mean:extended;
begin
  result:=r*(1-p)/p;
end;

function tNegative_binomial_distribution.variance:extended;
begin
  result:=mean/p;
end;

{-------------------------------------------------------------------}
{  tGeometric_distribution                                          }
{-------------------------------------------------------------------}

constructor tGeometric_distribution.create(newp:extended);
begin
  Inherited create(1,newp);
end;

{-------------------------------------------------------------------}
{  Initialization                                                   }
{-------------------------------------------------------------------}

begin
  randomize;
  seed1:=round(random*maxint)+1;
  seed2:=round(random*maxint)+1;
end.
