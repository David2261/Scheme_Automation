unit tau2;  {23.05.06}
interface
uses sysutils,dialogs;
const
  maxn = 10; {max stepen' PF}
  maxtabl = 100; {max chislo strok v tabl}
  maxmatr = 5; {K-vo tochek dlya D-razbieniya 30}
  PFAccuracy = 0.00000001; {Tochnost vichisleniy}
  DblAccuracy = 1e-30;
  MaxMemoPoints = 254; {Dlya ob'ekta BlockT}

type
  PFVector = array[0..maxn] of real;
  PFByteVector = array[0..maxn] of byte;
  TablPK = record
      n : byte;
      t, y : array[0..maxtabl] of real
    end;
  TablAFH = record
      n : byte;
      w, Re, Im, A, fi : array[0..maxtabl] of real
    end;
  TablD = record
      n : byte;
      w, k0, k1, k2 : array[0..maxtabl] of real
    end;

  {Transfer function:}
  TransferFunction = record
      a, b : PFVector;
      n, m : byte;
      k : real;
      tau : real;
    end;    
  PFArray = array[1..maxn] of TransferFunction;

  {Complex variable: (09.06.2004)}
  complex = object
      Re,Im : real;
      procedure S(mdl,fi : real);  {set by mod and arg}
      procedure E(a : complex);   {=}
      procedure P(a,b : complex); {+}
      procedure M(a,b : complex); {-}
      procedure U(a,b : complex); {*}
      procedure D(a,b : complex); {/}
      procedure K(a : real);      {multiply to real}
      procedure Exponent(a : complex);
      procedure Poly(a : PFVector; n : byte; x : complex);
      function Modul : real;  {modul}
      function F : real;     {argument}
      function stroka : string;
      end;

  PFComplexVector = array[0..maxn] of complex;
  dmatr = array[0..maxmatr,0..maxmatr] of real; {Promejutochnyy tip dlya D-razbieniya}
  
  TypePPK = record tr, est, yust, A1, A3, psi, sigma : real end; {Pryamye pokazateli}
  TypeKPK = record m, n : real; ust : boolean end; {Kornevye pokazateli. ust = true esli vse korni levye.}
  TypeFPK = record M, ME, wr, we, dA, dfi : real end; {Chastotnye pokazateli}

  {Blochnye PF:}
  block2 = object
      x1, y1, x2, y2,   {predyduschie znacheniya}
      c10, c11, c12, c21, c22 : real;
      procedure init(k, a0, a1, a2, b0, b1, b2, x10, y10, x20, y20 : real);
      procedure change(k, a0, a1, a2, b0, b1, b2 : real);
      function run(inp : real) : real;
    end;
  blockN = object  {PF N-stepeni}
      f1,f2 : PFVector;   {Koeffs for run}
      a,b,xt,yt : PFVector; {Input parameters}
      k : real; {Koeff}
      procedure zeros;
      procedure init;
      function run(x : real) : real;
    end;
  {Block overdue:}
  arr = array[1..MaxMemoPoints] of real;
  parr = ^arr;
(*  blockT = object
      t , {tau}
      cur_outp : real;
      n, dn, cur_t, cur_dt : byte;
      memory : parr;
      procedure init(tau, start : real);
      procedure change(tau, start : real);
      function run(inp : real) : real;
      procedure kill;
    end;*)
  blockT = object
      t , {tau}
      cur_outp : real;
      n, dn, cur_t, cur_dt : byte;
      memory : arr;
      procedure init(tau, start : real);
      procedure change(tau, start : real);
      function run(inp : real) : real;
    end;
  matr = array[0..maxn,0..maxn] of real;

var
{  CurrentPF : TransferFunction; {Rassmatrivaemaya PF objekta}
{  CurrentY0 : PFVector; {Nach. usloviya}
{  CurrentY0n : byte;
{  CurrentTablPK : TablPK;}
  CurrentdX : real; {Amplituda vhodnogo signala}
  SPF : PFArray; {Simou}
  SPFn : byte;   {Simou}
  Sm,Ss : PFVector; {Simou}
  CurrentL : byte; {Chislo korney (Laplas)}
  TimeDiscr : real = 0.1;  {Diskretnost po vremeni}



  
function Otklonenie(h1, h2 : TablPK) : real; {Raschet integralnogo otkloneniya}
procedure GetParmsFromStr(s : string; var x : PFVector; var nx : byte); {Chtenie stroki chisel}
procedure PolyD(a : PFVector; n : byte; b : PFVector; m : byte; var c : PFVector; var l : byte); {Delenie polinomov}
procedure PolyU(a : PFVector; n : byte; b : PFVector; m : byte; var c : PFVector; var l : byte); {Umnojenie polinomov}
procedure Roots(a : PFVector; m : byte; var x : PFComplexVector; var n : byte); {Korni}
{Laplas:}
procedure FindYs(pf : TransferFunction; vhod : real; y0 : PFVector; y0n : byte; var ys : TransferFunction); {1}
procedure LaplasM(pf : TransferFunction;
    var x, M : PFComplexVector; var l : byte; var compl, krat, nomer : PFByteVector); {2}
procedure LaplasY(x, M : PFComplexVector; l : byte; compl, krat, nomer : PFByteVector;
    tau, t0, dt : real; nt : byte; var PK : TablPK); {3}
function LaplasS(x, M : PFComplexVector; l : byte; compl, krat, nomer : PFByteVector; tau : real) : string; {4}
procedure Laplas(pf : TransferFunction; vhod : real; y0 : PFVector; y0n : byte; t0,dt : real; nt : byte;
    var PK : TablPK; var s : string); {0}
{AFH:}
procedure afhwm(ws : TransferFunction; m : real; w0, dw : real; nw : integer; typ : byte;
    var tabl : TablAFH); {Raschet AFH po zad. PF i m ili n.}
procedure afhasr(wr,wob : TransferFunction; m : real; w0, dw : real; nw : integer; typ : byte;
    var tablWraz, tablD, tablFz, tablFe, tablFv : TablAFH); {Raschet AFH po zad. PF regulatorai objecta i m ili n.}
{D-razbienie:}
procedure dmn(ws : TransferFunction; m : real; typ : byte; alfa : real; w0, dw : real; nw : byte;
    var tabl : TablD); {Krivaya D-razb dlya: typ=0 - raschet dlya m, typ=1 - dlya nu}
procedure dmme(ws : TransferFunction; m : real; typ : byte; alfa : real; w0, dw : real; nw : byte;
    var tabl : TablD); {Krivaya D-razb dlya: typ=0 - raschet dlya M, typ=1 - dlya ME}
{Simou:}
procedure Simou(h : TablPK; inp : real; var PF : PFArray; var PFn : byte; var m,s : PFVector);
{PK:}
procedure FindPPK(PK : TablPK; x : real; var PPK : TypePPK; var err : byte); {Poisk pryamyh pokazateley. err=1 - process ne ustanovilsya, err=2 - shag po vremini <=0}
procedure FindKPK(pf : TransferFunction; var KPK : TypeKPK; var err : byte); {Poisk kornevyh pokazateley. err=1 - net korney}
procedure FindFPK(tablWraz, tablFz, tablFe : TablAFH; var FPK : TypeFPK; var err : byte); {Poisk chastotnyh pokazateley. tablWraz, tablFz, tablFe - AFH Wraz, Fz i Fe}
{Perehodnye krivye:}
procedure ASRProcesses(wr,wob : TransferFunction; y0 : PFVector; y0n : byte; x : real; t0, dt : real; nt : byte;
    var yt, et, ut : TablPK); {Raschet perehodnyh processov v ASR: wr, wob - PF regulyatora i OU, y0 - vector nach usloviy dlya OU}



implementation

{BEGIN - Service Functions ------------------------------------------------------}
function pow(k:integer;a:real):real; {Used in Simou}
var p:real;
begin
 if a=0 then p:=0
 else
 if k=0 then p:=1
 else
  begin
   p:=exp(k*ln(abs(a)));
   if a<0 then
    if (k mod 2)<>0 then p:=-p;
  end;
 pow:=p;
end{pow};

function Otklonenie(h1, h2 : TablPK) : real;
{Raschet integralnogo otkloneniya h2(t) ot h1(t)}
{Kolichestvo tochek v h1 i h2 doljno byt' odinakovym = n}
var i,n : integer; res : real;
begin
  res:=0; if h1.n>h2.n then n:=h2.n else n:=h1.n;
  for i:=0 to n do res:=res+(h1.y[i]-h2.y[i])*(h1.y[i]-h2.y[i]);
  Otklonenie:=res
end;

(*procedure GetParmsFromStr(s : string; var x : PFVector; var nx : byte);
{Chitaet stroku i videlyaet iz nee chisla. Bukvy i lishnie '.' ignoriruet.}
{nx - chislo parametrov}
var i,n : byte; s1 : string;  cod : integer; decsep : boolean;
begin
   {DecimalSeparator:='.';}
    if s='' then begin nx:=0; exit end; 
    n:=length(s); decsep:=false; nx:=0; s1:=''; {i - nomer $, nx - nomer peremennoy}
    for i:=0 to n do
      	case s[i] of
      	'0'..'9': s1:=s1+s[i];
      	',':  if not decsep then  {DecimalSeparator}
      		     begin s1:=s1+'.'; decsep:=true end;
      	'.':  if not decsep then  {DecimalSeparator}
      		     begin s1:=s1+'.'; decsep:=true end;
        ' ': if s1<>'' then begin val(s1,x[nx],cod); nx:=nx+1; s1:=''; decsep:=false end;
        '-': if s1='' then s1:=s1+'-'
        end;
    if (s1<>'') and (s1<>'-') then begin val(s1,x[nx],cod); nx:=nx+1 end
end;             *)

procedure GetParmsFromStr(s : string; var x : PFVector; var nx : byte);
{Chitaet stroku i videlyaet iz nee chisla. Bukvy i lishnie '.' ignoriruet.}
{nx - chislo parametrov}
var i,n : byte; s1 : string;  cod : integer; decsep : boolean;
begin
   {DecimalSeparator:='.';}
    if s='' then begin nx:=0; exit end; 
    n:=length(s); decsep:=false; nx:=0; s1:=''; {i - nomer $, nx - nomer peremennoy}
    for i:=0 to n do
      	case s[i] of
      	'0'..'9': s1:=s1+s[i];
      	',':  if not decsep then  {DecimalSeparator}
      		     begin s1:=s1+DecimalSeparator; decsep:=true end;
      	'.':  if not decsep then  {DecimalSeparator}
      		     begin s1:=s1+DecimalSeparator; decsep:=true end;
        ' ': if s1<>'' then begin x[nx]:=strtofloat(s1); nx:=nx+1; s1:=''; decsep:=false end;
        '-': if s1='' then s1:=s1+'-'
        end;
    if (s1<>'') and (s1<>'-') then begin x[nx]:=strtofloat(s1); nx:=nx+1 end
end;

{END - Service Functions ------------------------------------------------------}






{BEGIN - Complex variable ------------------------------------------------------}
procedure complex.S;
  begin Re:=mdl*cos(fi/180*pi); Im:=mdl*sin(fi/180*pi) end;
procedure complex.E;  begin Re:=a.Re; Im:=a.Im end;
procedure complex.P;  begin Re:=a.Re+b.Re; Im:=a.Im+b.Im end;
procedure complex.M;  begin Re:=a.Re-b.Re; Im:=a.Im-b.Im end;
procedure complex.U;
  begin Re:=a.Re*b.Re-a.Im*b.Im; Im:=a.Re*b.Im+a.Im*b.Re end;
procedure complex.D;
  var r:real;
  begin r:=b.Re*b.Re+b.Im*b.Im;
    if r=0 then Re:=0 else Re:=(a.Re*b.Re+a.Im*b.Im)/r;
    if r=0 then Im:=0 else Im:=(a.Im*b.Re-a.Re*b.Im)/r end;
procedure complex.K(a : real);      {multiply to real}
  begin Re:=a*Re; Im:=a*Im end;
procedure complex.Exponent(a : complex);
  var ex : real;
  begin ex:=exp(a.Re); Re:=ex*cos(a.Im); Im:=ex*sin(a.Im) end;
procedure complex.Poly(a : PFVector; n : byte; x : complex);
  var sum,q,xx : complex; i : byte;
  begin
    sum.re:=a[0]; sum.im:=0;
    xx.e(x);
    for i:=1 to n do
      begin q.e(xx); q.k(a[i]); sum.p(sum,q); xx.u(xx,x) end;
    Re:=sum.re; Im:=sum.im
  end;
function complex.modul : real;
  var r:real; begin r:=sqrt(Re*Re+Im*Im); modul:=r end;
function complex.F : real;
  var fi:real;
  begin if Re=0 then if Im>=0 then fi:=90 else fi:=-90
    else fi:=arctan(Im/Re)/pi*180;
    if Re<0 then fi:=fi+180; f:=fi end;
function complex.stroka : string;
  var st,t : string; r:real;
  begin str(Re:8:4,t); st:='('+t+')+j*('; str(Im:8:4,t); st:=st+t;
    r:=modul; str(r:8:4,t); st:=st+') = ('+t+')*e^j(';
    r:=f; str(r:8:4,t); st:=st+t+')'; stroka:=st end;

function IsEqual(x1, x2 : complex) : boolean;
{=true, if x1=x2}
  var res : boolean;
  begin
    res:=false;
    if (abs(x1.re-x2.re)<=PFAccuracy) and (abs(x1.im-x2.im)<=PFAccuracy) then res:=true;
    ISEqual:=res
  end;
{END - Comlex variable ------------------------------------------------------}






{BEGIN - Polynoms ------------------------------------------------------}
procedure PolyD(a : PFVector; n : byte; b : PFVector; m : byte; var c : PFVector; var l : byte);
{Delenie polinoma A na polinom B. Resultat: polinom C.}
var i,j : byte;
begin
  {ubiraem nuli iz B}
  j:=0; for i:=0 to m do if abs(b[i])>PFAccuracy then j:=i;
  m:=j;
  c[0]:=0;
  if n<m then begin l:=0; exit end;
  l:=n-m;
  for i:=l downto 0 do
    begin
      c[i]:=a[n]/b[m];
      for j:=0 to m do a[n-j]:=a[n-j]-c[i]*b[m-j];
      n:=n-1
    end
end;

procedure PolyU(a : PFVector; n : byte; b : PFVector; m : byte; var c : PFVector; var l : byte);
{Umnojenie polinoma A na polinom B. Resultat: polinom C.}
var i,j : byte;
begin
  {ubiraem nuli iz B}
  j:=0; for i:=0 to m do if abs(b[i])>PFAccuracy then j:=i;
  m:=j;
  l:=m+n;
  for i:=0 to l do c[i]:=0;
  for i:=0 to m do
      for j:=0 to n do
          c[i+j]:=c[i+j]+a[j]*b[i]
end;
{END - Polynoms ------------------------------------------------------}







{BEGIN - Block PFs ------------------------------------------------------}
{Block 2: W(s) = K*(b2*s2 + b1*s + b0) / (a2*s + a1*s + a0) }
procedure block2.init(k, a0, a1, a2, b0, b1, b2, x10, y10, x20, y20 : real);
  var p : real;
  begin
    p:=a0*TimeDiscr*TimeDiscr+a1*TimeDiscr+a2;
    c10:=k*(b0*TimeDiscr*TimeDiscr+b1*TimeDiscr+b2)/p;
    c11:=-k*(b1*TimeDiscr+2*b2)/p;
    c12:=k*b2/p;
    c21:=(a1*TimeDiscr+2*a2)/p;
    c22:=-a2/p;
    x1:=x10; y1:=y10; x2:=x20; y2:=y20
  end;
procedure block2.change(k, a0, a1, a2, b0, b1, b2 : real);
  var p : real;
  begin
    p:=a0*TimeDiscr*TimeDiscr+a1*TimeDiscr+a2;
    c10:=k*(b0*TimeDiscr*TimeDiscr+b1*TimeDiscr+b2)/p;
    c11:=-k*(b1*TimeDiscr+2*b2)/p;
    c12:=k*b2/p;
    c21:=(a1*TimeDiscr+2*a2)/p;
    c22:=-a2/p;
  end;
function block2.run(inp : real) : real;
  var r : real;
  begin
    r:=c10*inp+c11*x1+c12*x2+c21*y1+c22*y2;
    x2:=x1; y2:=y1; x1:=inp; y1:=r; run:=r
  end;

{Block N: W(s) = K*(bn*sn+...+b2*s2 + b1*s + b0) / (an*sn+...+a2*s2 + a1*s + a0) }
procedure blockN.init;
  var
    c,d1,d2 : matr;
    T : PFVector; {Stepeni TimeDiscr}
    i,j : byte;
    r,sum : real;
  begin
    if TimeDiscr=0 then exit;
    {Calc C:}
    for i:=0 to maxn do c[0,i]:=0;
    for i:=0 to maxn do c[i,0]:=1;
    for i:=1 to maxn do
        for j:=1 to maxn do
            c[i,j]:=c[i-1,j]-c[i-1,j-1];
    {Calc D1,D2:}
    r:=1; sum:=a[0]; T[0]:=1;
    for i:=1 to maxn do begin r:=r*TimeDiscr; T[i]:=r; sum:=sum+a[i]/r end;
    if sum=0 then exit;
    for i:=0 to maxn do
        for j:=0 to maxn do
          begin
            d1[i,j]:=k*b[i]/T[i]*c[i,j]/sum;
            d2[i,j]:=a[i]/T[i]*c[i,j]/sum
          end;
    {Calc F1, F2:}
    for j:=0 to maxn do
      begin
        f1[j]:=0; f2[j]:=0;
        for i:=0 to maxn do
          begin f1[j]:=f1[j]+d1[i,j]; f2[j]:=f2[j]+d2[i,j] end
      end
  end;

function blockN.run(x : real) : real;
  var i : byte; y : real;
  begin
    for i:=maxn downto 1 do
      begin xt[i]:=xt[i-1]; yt[i]:=yt[i-1] end;
    y:=f1[0]*x;
    for i:=1 to maxn do y:=y+f1[i]*xt[i]-f2[i]*yt[i];
    xt[0]:=x; yt[0]:=y;
    run:=y
  end;

procedure blockN.zeros; {obnulenie nach. znacheniy}
  var i:byte;
  begin
    for i:=0 to maxn do
      begin a[i]:=0; b[i]:=0; xt[i]:=0; yt[i]:=0; f1[i]:=0; f2[i]:=0 end;
    k:=1
  end;

{Block T: W(s) = exp(-T*s)  }
procedure blockT.init(tau, start : real);
  var m,i : integer;
  begin
    t:=tau; n:=0; if t<=0 then exit;
    m:=round(t/TimeDiscr); cur_t:=0; cur_dt:=0; cur_outp:=0;
    if m<MaxMemoPoints then
      begin dn:=round(t/TimeDiscr/MaxMemoPoints+0.51); n:=round(t/TimeDiscr/dn) end
      else begin dn:=1; n:=m end;
    for i:=1 to n do memory[i]:=start;
  end;
procedure blockT.change(tau, start : real);
  var m,i : integer;
  begin
{MessageDlg('1',mtInformation,[mbOk],0);}
{    for i:=1 to n do memory^[i]:=start;}
    {if SEqu(t,tau) then exit;}
    if tau<=0 then exit;
    t:=tau; n:=0;
    {if memory<>nil then Dispose(memory);
    New(memory);}
    m:=round(t/TimeDiscr); cur_t:=0; cur_dt:=0; cur_outp:=start;
    if m<MaxMemoPoints then
      begin dn:=round(t/TimeDiscr/MaxMemoPoints+0.51); n:=round(t/TimeDiscr/dn) end
      else begin dn:=1; n:=m end;
    for i:=1 to n do memory[i]:=start;
  end;
function blockT.run(inp : real) : real;
  begin
    if n=0 then cur_outp:=inp
      else begin
        cur_dt:=cur_dt+1;
        if cur_dt=dn then {pora peredavat/pisat znachenie}
          begin
            cur_dt:=0; cur_t:=cur_t+1;  if cur_t>n then cur_t:=1;
            cur_outp:=memory[cur_t];
            memory[cur_t]:=inp
          end
      end;
    run:=cur_outp
  end;
(*procedure blockT.kill;  begin Dispose(memory) end;*)
{END - Block PFs ------------------------------------------------------}








{BEGIN - Roots ------------------------------------------------------}
procedure Roots(a : PFVector; m : byte; var x : PFComplexVector; var n : byte);
{A - polinom, m - razmer A, X - korni, n - chislo korney}
var i,l : byte; a10,b10,a20,b20,x1,x2,y,d,a2 : real; p : boolean; b : PFVector;
label lbl;

function f(x1,x2:real):real; {Function for look}
{minimiziruemaya funkciya}
  var xx,c : complex;
  begin
    xx.re:=x1; xx.im:=x2; c.poly(a,m,xx);
    f:=c.re*c.re+c.im*c.im
  end;

procedure FindZero(a10,b10,a20,b20 : real);
var
  a1,b1,a2,b2,x11,x12,x21,x22,epsn,y,x1_itog,x2_itog : real;
  i,j : byte;
  fy : array[1..4] of real;
begin
  a1:=a10; b1:=b10; a2:=a20; b2:=b20;
  repeat
    x11:=0.6*a1+0.4*b1;
    x12:=0.4*a1+0.6*b1;
    x21:=0.6*a2+0.4*b2;
    x22:=0.4*a2+0.6*b2;
    fy[1]:=f(x12,x22);
    fy[2]:=f(x11,x22);
    fy[3]:=f(x11,x21);
    fy[4]:=f(x12,x21);
    j:=1;
    for i:=1 to 4 do if fy[i]>fy[j] then j:=i;
    case j of
    1: begin b1:=x12; b2:=x22 end;
    2: begin a1:=x11; b2:=x22 end;
    3: begin a1:=x11; a2:=x21 end;
    4: begin b1:=x12; a2:=x21 end;
    end;
    epsn:=(b1-a1)/2;
  until epsn<=PFAccuracy;
  x1_itog:=(a1+b1)/2; x2_itog:=(a2+b2)/2;
  {Analiz resheniya:}
  y:=f(x1_itog,x2_itog);
  if y<0.6*PFAccuracy then
    begin x1:=x1_itog; x2:=x2_itog; p:=true end
    else begin {eto esche ne koren'}
      if p then exit;
      FindZero(x1_itog,b10,x2_itog,b20);
      FindZero(a10,x1_itog,x2_itog,b20);
      FindZero(a10,x1_itog,a20,x2_itog);
      FindZero(x1_itog,b10,a20,x2_itog);
    end
end;

begin
  n:=0;
  {vydelenie nuley:}
  p:=true;
  for i:=0 to m-1 do
      if (p) and (abs(a[i])<PFAccuracy) then {nul}
        begin n:=n+1; x[n].re:=0; x[n].im:=0 end
        else p:=false;
  if n>0 then for i:=0 to m-n do a[i]:=a[i+n];
  if n>m then n:=0 else m:=m-n;
  {Opredelenie diapazona poiska:}
  y:=a[0];
  for i:=0 to m do if y<abs(a[i]) then y:=abs(a[i]);
  a10:=-2*y-10*PFAccuracy; b10:=2*y+10*PFAccuracy; a20:=0; b20:=2*y+10*PFAccuracy;
lbl:
    if m=1 then begin n:=n+1; x[n].re:=-a[0]/a[1]; x[n].im:=0 end;
    if m=2 then
      begin
        d:=a[1]*a[1]-4*a[0]*a[2]; {diskriminant}
        a2:=2*a[2];
        if d>0 then
          begin
            d:=sqrt(d);
            n:=n+1; x[n].re:=(-a[1]+d)/a2; x[n].im:=0;
            n:=n+1; x[n].re:=(-a[1]-d)/a2; x[n].im:=0
          end
          else begin
            d:=sqrt(-d);
            n:=n+1; x[n].re:=-a[1]/a2; x[n].im:=d/a2;
            n:=n+1; x[n].re:=x[n-1].re; x[n].im:=-x[n-1].im
          end
      end;
    if m>2 then
    begin {Poisk korney:}
      p:=false;
      FindZero(a10,b10,a20,b20);
      y:=f(x1,x2);
      if (y<0.6*PFAccuracy) then {It's a root}
        begin
          if n=maxn then exit;
          n:=n+1; x[n].re:=x1; x[n].im:=x2;
          if n=maxn then exit;
          if abs(x2)>PFAccuracy then {It's a complex}
            begin
              n:=n+1; x[n].re:=x1; x[n].im:=-x2;
              b[2]:=1; b[1]:=-2*x1; b[0]:=x1*x1+x2*x2; l:=2
            end
            else begin
              b[1]:=1; b[0]:=-x1; l:=1
            end;
          PolyD(a,m,b,l,a,m);
          goto lbl
        end
    end;
end;
{END - Roots ------------------------------------------------------}






{BEGIN - Laplas ------------------------------------------------------}
procedure LaplasM(pf : TransferFunction;
var x, M : PFComplexVector; var l : byte; var compl, krat, nomer : PFByteVector);
{Poisk k-tov Mi. k,a,n,b,m,tau - parametry PF; x- korni, M- k-ty, l - chislo korney}
  var i,j,kratnomer : byte; q : complex;

  procedure FindM(i : byte); {Raschet Mi. i - nomer kornya}
  {!!! - Ne realizovan raschet dlya kratnyh korney!}
    var y1,y2 : complex; c : PFVector; ii : byte;
    begin
      {Chislitel:}
      y1.poly(pf.b,pf.n,x[i]); y1.k(pf.k);
      {Znamenatel:}
      for ii:=0 to pf.n-1 do c[ii]:=(ii+1)*pf.a[ii+1];
      y2.poly(c,pf.n-1,x[i]);
      {Koefficient:}
      m[i].d(y1,y2);
    end;

  begin
    {Korni:}
    Roots(pf.a,pf.n,x,l);
    if l=0 then exit;
    {Analiz korney:}
    for i:=1 to maxn do
      begin M[i].re:=0; M[i].im:=0; compl[i]:=0; krat[i]:=0; nomer[i]:=0 end;
    for i:=1 to l do
      begin
        if abs(x[i].im)>PFAccuracy then {It's a complex}
          begin
            compl[i]:=1;
            if nomer[i]=0 then
              begin
                nomer[i]:=1;
                q.re:=x[i].re; q.im:=-x[i].im;
                {Raschet Mi:}
                FindM(i);
                {Poisk parnogo:}
                for j:=i+1 to l do
                    if IsEqual(x[j],q) then
                      begin nomer[j]:=2; M[j].re:=M[i].re; M[j].im:=-M[i].im end
              end
          end
          else begin {It's a real}
            {Proverka na kratnost':}
            kratnomer:=1;
            if krat[i]=0 then
              begin
                for j:=i+1 to l do
                    if IsEqual(x[i],x[j]) then
                      begin
                        kratnomer:=kratnomer+1;
                        krat[i]:=1; krat[j]:=1;
                        nomer[i]:=1; nomer[j]:=kratnomer;
                        {Raschet Mi:}
                        FindM(i); FindM(j)
                      end
              end;
            if krat[i]=0 then FindM(i)
          end
      end;
  end;

procedure LaplasY(x, M : PFComplexVector; l : byte; compl, krat, nomer : PFByteVector;
tau, t0, dt : real; nt : byte; var PK : TablPK);
{Raschet perehodnoy krivoy po x[i], M[i]. Vyhod: tablica PK}
  var i,j : byte; t : real;
  begin
    if nt=0 then exit;
    t:=t0; PK.n:=nt;
    for i:=1 to nt do
      begin
        PK.t[i]:=t; PK.y[i]:=0;
        if t>=tau then
          begin
            for j:=1 to l do
              begin
                if (compl[j]=0) and (krat[j]=0) then {Koren deystv. ne kratn.}
                    PK.y[i]:=PK.y[i]+M[j].re*exp(x[j].re*(t-tau));
                if (compl[j]=1) and (nomer[j]=1) then {Koren compl. #1}
                    PK.y[i]:=PK.y[i]+2*exp(x[j].re*(t-tau))*
                    (M[j].re*cos(x[j].im*(t-tau))-M[j].im*sin(x[j].im*(t-tau)))
              end
          end;
        t:=t+dt
      end;
  end;

function LaplasS(x, M : PFComplexVector; l : byte; compl, krat, nomer : PFByteVector; tau : real) : string;
{Predstavlenie perehodnoy harakteristiki}
  var j : byte; s,st,p,p1,p2,p3,p4 : string;
  begin
    s:='y(t) =';
    str(tau:4:0,p);
    if tau<PFAccuracy then st:='t' else st:='(t-'+p+')';
    for j:=1 to l do
      begin
        if (compl[j]=0) and (krat[j]=0)  then {Koren deystv. ne kratn.}
          begin
            if abs(x[j].re)<PFAccuracy then {nulevoy koren'}
              begin
                str(M[j].re:8:4,p);
                if (j<>1) and (M[j].re>0) then s:=s+' +';
                s:=s+' '+p
              end
              else begin
                str(M[j].re:8:4,p1); str(x[j].re:8:4,p2);
                if (j<>1) and (M[j].re>0) then s:=s+' +';
                s:=s+' '+p1+'*exp('+p2+'*'+st+')'
              end
          end;
        if (compl[j]=1) and (nomer[j]=1) then {Koren compl. #1}
          begin
            if j<>1 then s:=s+' +';
            if abs(x[j].re)<PFAccuracy then {Re=0}
              begin
                str(M[j].re:8:4,p1); str(M[j].im:8:4,p2); str(x[j].im:8:4,p);
                s:=s+' 2*('+p1+'*cos('+p+'*'+st+') ';
                if M[j].im<0 then s:=s+'+' else s:=s+'-'; {znaki}
                s:=s+' '+p2+'*sin('+p+'*'+st+'))'
              end
              else begin
                str(M[j].re:8:4,p1); str(M[j].im:8:4,p2); str(x[j].re:8:4,p3); str(x[j].im:8:4,p4);
                s:=s+' 2*exp('+p3+'*'+st+')*('+p1+'*cos('+p4+'*'+st+') ';
                if M[j].im<0 then s:=s+'+' else s:=s+'-'; {znaki}
                s:=s+' '+p2+'*sin('+p4+'*'+st+'))'
              end
          end;
      end;
    LaplasS:=s
  end;

procedure FindYs(pf : TransferFunction; vhod : real; y0 : PFVector; y0n : byte; var ys : TransferFunction);
{Formiruet Y(s) po PF s uchetom nachal'nyh usloviy (X(s) = vhod/s)}
  var i, j : byte;
  begin
    ys.k:=1; ys.tau:=pf.tau; ys.n:=pf.n+1; ys.a[0]:=0;
    for i:=0 to pf.n do ys.a[i+1]:=pf.a[i];
    if y0n=0 then {n.n.u.}
      begin
        ys.m:=pf.m;
        for i:=0 to pf.m do ys.b[i]:=pf.k*pf.b[i]*vhod
      end
      else begin
      	ys.m:=pf.n;
        if pf.m<maxn then for i:=pf.m+1 to maxn do pf.b[i]:=0;
      	for i:=0 to pf.n do 
      	  begin
      	    ys.b[i]:=pf.k*pf.b[i]*vhod;
      	    for j:=i to pf.n do ys.b[i]:=ys.b[i]+pf.a[j]*y0[j-i]
      	  end
      end
  end;

procedure Laplas(pf : TransferFunction; vhod : real; y0 : PFVector; y0n : byte; t0,dt : real; nt : byte;
var PK : TablPK; var s : string);
{Raschet perehodnoy harakteristiki po PF i nach. usloviyam}
  var x, M : PFComplexVector; l : byte; compl, krat, nomer : PFByteVector; ys : TransferFunction;
  begin
    FindYs(pf,vhod,y0,y0n,ys);
    LaplasM(ys,x,M,l,compl,krat,nomer);
    LaplasY(x,M,l,compl,krat,nomer,pf.tau,t0,dt,nt,PK);
    s:=LaplasS(x,M,l,compl,krat,nomer,pf.tau);
  end;

{END - Laplas ------------------------------------------------------}





{BEGIN - AFH ------------------------------------------------------}
procedure afhwm(ws : TransferFunction; m : real; w0, dw : real; nw : integer; typ : byte;
var tabl : TablAFH);
{Raschet AFH po zad. PF i m ili n.}
{typ=0 - raschet dlya m; typ=1 - raschet dlya n}
{Resultat v tabl}
  var
    i,j : integer;
    w : real;
    s,p1,p2,a,b,e : complex;
  begin
    w:=w0;
    for j:=1 to nw do
      begin
        if typ=0 then s.re:=-m*w else s.re:=m;
        s.im:=w;
        {B(s):}
        b.re:=ws.b[0]; b.im:=0; p1.e(s);
        for i:=1 to ws.m do
          begin
            p2.e(p1); p2.k(ws.b[i]); b.p(b,p2); p1.u(p1,s)
          end;
        b.k(ws.k);
        {A(s):}
        a.re:=ws.a[0]; a.im:=0; p1.e(s);
        for i:=1 to ws.n do
          begin
            p2.e(p1); p2.k(ws.a[i]); a.p(a,p2); p1.u(p1,s)
          end;
        {e(s):}
        p1.e(s); p1.k(-ws.tau); e.exponent(p1);
        p1.u(b,e); p2.d(p1,a);
        {Tabl:}
        tabl.w[j]:=w; tabl.re[j]:=p2.re; tabl.im[j]:=p2.im;
        tabl.a[j]:=p2.modul; tabl.fi[j]:=p2.f;
        w:=w+dw
      end;
    tabl.n:=nw
  end;

procedure afhasr(wr,wob : TransferFunction; m : real; w0, dw : real; nw : integer; typ : byte;
var tablWraz, tablD, tablFz, tablFe, tablFv : TablAFH);
{Raschet AFH po zad. PF regulatorai objecta i m ili n.}
{pftyp=0 - Wraz(s), pftyp=1 - D(s), pftyp=2 - Fz(s), pftyp=3 - Fe(s), pftyp=4 - Fv(s)}
{typ=0 - raschet dlya m; typ=1 - raschet dlya n}
{Resultat v tabl}
  var
    i,j : integer;
    w : real;
    s,p1,p2,ar,br,aob,bob,b,e,pr,pob,praz,res : complex;
  begin
    w:=w0;
    for j:=1 to nw do
      begin
        if typ=0 then s.re:=-m*w else s.re:=m;
        s.im:=w;

        {Regulator Wr(s):}
        {Br(s):}
        b.re:=wr.b[0]; b.im:=0; p1.e(s);
        if wr.m>0 then for i:=1 to wr.m do
          begin
            p2.e(p1); p2.k(wr.b[i]); b.p(b,p2); p1.u(p1,s)
          end;
        b.k(wr.k);
        {Ar(s):}
        ar.re:=wr.a[0]; ar.im:=0; p1.e(s);
        if wr.n>0 then for i:=1 to wr.n do
          begin
            p2.e(p1); p2.k(wr.a[i]); ar.p(ar,p2); p1.u(p1,s)
          end;
        {e(s):}
        p1.e(s); p1.k(-wr.tau); e.exponent(p1);
        br.u(b,e); pr.d(br,ar);

        {Object Wob(s):}
        {Bob(s):}
        b.re:=wob.b[0]; b.im:=0; p1.e(s);
        if wob.m>0 then for i:=1 to wob.m do
          begin
            p2.e(p1); p2.k(wob.b[i]); b.p(b,p2); p1.u(p1,s)
          end;
        b.k(wob.k);
        {Aob(s):}
        aob.re:=wob.a[0]; aob.im:=0; p1.e(s);
        if wob.n>0 then for i:=1 to wob.n do
          begin
            p2.e(p1); p2.k(wob.a[i]); aob.p(aob,p2); p1.u(p1,s)
          end;
        {e(s):}
        p1.e(s); p1.k(-wob.tau); e.exponent(p1);
        bob.u(b,e); pob.d(bob,aob);

(*        case pftyp of
        0: begin {Wraz}
             res.u(pr,pob)
           end;
        1: begin {D(s)}
             p1.u(br,bob); p2.u(ar,aob);
             res.p(p1,p2)
           end;
        2: begin {Fz(s)}
             p1.u(pr,pob); p2.e(p1); p2.re:=p2.re+1;
             res.d(p1,p2)
           end;
        3: begin {Fe(s)}
             p1.u(pr,pob); p2.re:=1; p2.im:=0; p1.re:=p1.re+1;
             res.d(p2,p1)
           end;
        4: begin {Fv(s)}
             p1.u(pr,pob); p1.re:=p1.re+1;
             res.d(pob,p1)
           end;
        end;

        {Tabl:}
        tabl.w[j]:=w; tabl.re[j]:=res.re; tabl.im[j]:=res.im;
        tabl.a[j]:=res.modul; tabl.fi:=res.f[j]; *)
        
        {Wraz}
        res.u(pr,pob);
        tablWraz.w[j]:=w; tablWraz.re[j]:=res.re; tablWraz.im[j]:=res.im;
        tablWraz.a[j]:=res.modul; tablWraz.fi[j]:=res.f;
        {D(s)}
        p1.u(br,bob); p2.u(ar,aob);res.p(p1,p2);
        tablD.w[j]:=w; tablD.re[j]:=res.re; tablD.im[j]:=res.im;
        tablD.a[j]:=res.modul; tablD.fi[j]:=res.f;
        {Fz(s)}
        p1.u(pr,pob); p2.e(p1); p2.re:=p2.re+1; res.d(p1,p2);
        tablFz.w[j]:=w; tablFz.re[j]:=res.re; tablFz.im[j]:=res.im;
        tablFz.a[j]:=res.modul; tablFz.fi[j]:=res.f;
        {Fe(s)}
        p1.u(pr,pob); p2.re:=1; p2.im:=0; p1.re:=p1.re+1; res.d(p2,p1);
        tablFe.w[j]:=w; tablFe.re[j]:=res.re; tablFe.im[j]:=res.im;
        tablFe.a[j]:=res.modul; tablFe.fi[j]:=res.f;
        {Fv(s)}
        p1.u(pr,pob); p1.re:=p1.re+1; res.d(pob,p1);
        tablFv.w[j]:=w; tablFv.re[j]:=res.re; tablFv.im[j]:=res.im;
        tablFv.a[j]:=res.modul; tablFv.fi[j]:=res.f;

        w:=w+dw
      end;
    tablWraz.n:=nw; tablD.n:=nw; tablFz.n:=nw; tablFe.n:=nw; tablFv.n:=nw
  end;
{END - AFH ------------------------------------------------------}




(*
{BEGIN - D-razbienie po n i m ------------------------------------------------------}
procedure dmn(ws : TransferFunction; m : real; typ : byte; alfa : real; w0, dw : real; nw : byte;
var tabl : TablD);
{Krivaya D-razb dlya:}
{typ=0 - raschet dlya m, typ=1 - dlya nu}
  var
    i,j : byte;
    w,k0,k1,e,h : real;
    a,b,s,p1,p2,pb,pw,c1,c2,c3,c4 : complex;
  label Loop;

function func(k0,k1 : real) : real; {Function for look}
  var r,f1,f2:real;
  begin
    if k0=0 then r:=10000 else r:=k1/k0;
    f1:=c1.re*k0+c3.re*k1+c4.re*k1*r+c2.re;
    f2:=c1.im*k0+c3.im*k1+c4.im*k1*r+c2.im;
    func:=f1*f1+f2*f2
  end;

function FindZeroK0(a0,b0,eps : real) : real;
var x1,x2,f1,f2,a,b,epsn : real;
begin
  a:=a0; b:=b0;
  repeat                                      
    x1:=0.6*a+0.4*b;
    x2:=0.4*a+0.6*b;
    f1:=func(x1,k1);
    f2:=func(x2,k1);
    if f1<=f2 then b:=x2 else a:=x1;
    epsn:=(b-a)/2;
  until epsn<=eps;
  FindZeroK0:=(a+b)/2
end;

function FindZeroK1(a0,b0,eps : real) : real;
var x1,x2,f1,f2,a,b,epsn : real;
begin
  a:=a0; b:=b0;
  repeat
    x1:=0.6*a+0.4*b;
    x2:=0.4*a+0.6*b;
    f1:=func(k0,x1);
    f2:=func(k0,x2);
    if f1<=f2 then b:=x2 else a:=x1;
    epsn:=(b-a)/2;
  until epsn<=eps;
  FindZeroK1:=(a+b)/2
end;

{-----------------------MAIN -----------------}
  begin
    w:=w0;
    for j:=1 to nw do
      begin
        if typ=0 then s.re:=-m*w else s.re:=m;
        s.im:=w;

        {B(s):}
        b.re:=ws.b[0]; b.im:=0; p1.e(s);
        for i:=1 to ws.m do
          begin
            p2.e(p1); p2.k(ws.b[i]); b.p(b,p2); p1.u(p1,s)
          end;
        b.k(ws.k);

        {A(s):}
        a.re:=ws.a[0]; a.im:=0; p1.e(s);
        for i:=1 to ws.n do
          begin
            p2.e(p1); p2.k(ws.a[i]); a.p(a,p2); p1.u(p1,s)
          end;

        {e(s):}
        p1.e(s); p1.k(-ws.tau); p2.exponent(p1);
        pb.u(b,p2); pw.d(pb,a);

        {Koeffs:}
        c1.e(pb); c2.u(s,a); c3.u(s,c1); c4.u(s,c3); c4.k(alfa);

        {Solve:}
        h:=1000000; e:=0.0001; k0:=0; k1:=0;
        repeat
          k0:=FindZeroK0(k0-h,k0+h,e);
          k1:=FindZeroK1(k1-h,k1+h,e);
          h:=h/10;
        until func(k0,k1)<=e;
        if (w=0) and (ws.k<>0) then k1:=-1/ws.k;

        {Tabl:}
        tabl.w[j]:=w; tabl.k0[j]:=k0; tabl.k1[j]:=k1;
        if k0<>0 then tabl.k2[j]:=alfa*k1*k1/k0 else tabl.k2[j]:=0;
        w:=w+dw
      end;
    tabl.n:=nw
end;
{END - D-razbienie po n i m ------------------------------------------------------}
*)


{BEGIN - D-razbienie po n i m ------------------------------------------------------}
procedure dmn(ws : TransferFunction; m : real; typ : byte; alfa : real; w0, dw : real; nw : byte;
var tabl : TablD);
{Krivaya D-razb dlya:}
{typ=0 - raschet dlya m, typ=1 - dlya nu}
  var
    i,j : byte;
    w,k0,k1,e,h : real;
    a,b,s,p1,p2,pb,pw,c1,c2,c3,c4 : complex;
    d1,d2,d3,e1,e2,f0,f1,f2,discr,pp : double;
  label Loop;
  begin
    w:=w0;
    for j:=1 to nw do
      begin
{        if typ=0 then s.re:=-m*w else s.re:=m;
        s.im:=w;}
{        if abs(w)<DblAccuracy then begin s.Im:=DblAccuracy; w:=DblAccuracy end;}
        if abs(w)>PFAccuracy then
          begin {w<>0}
            if typ=0 then s.re:=-m*w else s.re:=m;
            s.im:=w
          end
          else
          begin {w=0}
              s.Im:=dw/5;
              s.Re:=0
          end;

        {B(s):}
        b.re:=ws.b[0]; b.im:=0; p1.e(s);
        for i:=1 to ws.m do
          begin
            p2.e(p1); p2.k(ws.b[i]); b.p(b,p2); p1.u(p1,s)
          end;
        b.k(ws.k);

        {A(s):}
        a.re:=ws.a[0]; a.im:=0; p1.e(s);
        for i:=1 to ws.n do
          begin
            p2.e(p1); p2.k(ws.a[i]); a.p(a,p2); p1.u(p1,s)
          end;

        {e(s):}
        p1.e(s); p1.k(-ws.tau); p2.exponent(p1); {p2 = exp(-tau*s)}
        b.u(b,p2); {a-znamenatel, b-chislitel}

        {Koeffs:}
        c1.u(s,a); c2.u(b,s); c3.e(b); c4.u(b,s); c4.u(c4,s); c4.k(alfa);
        if abs(c4.re)<DblAccuracy then begin if c4.re>0 then c4.re:=DblAccuracy else c4.re:=-DblAccuracy end;
        if abs(c4.im)<DblAccuracy then begin if c4.im>0 then c4.im:=DblAccuracy else c4.im:=-DblAccuracy end;
        d1:=c1.im-c1.re*c4.im/c4.re;
        d2:=c2.im-c2.re*c4.im/c4.re;
        d3:=c3.im-c3.re*c4.im/c4.re;
        if abs(d2)<DblAccuracy then begin if d2>0 then d2:=DblAccuracy else d2:=-DblAccuracy end;
        e1:=-d1/d2; e2:=-d3/d2;
        f2:=c2.re*e2+c3.re+c4.re*e2*e2;
        f1:=c2.re*e1+2*c4.re*e1*e2+c1.re;
        f0:=c4.re*e1*e1;
        if abs(f2)<DblAccuracy then begin if f2>0 then f2:=DblAccuracy else f2:=-DblAccuracy end;
        discr:=f1*f1-4*f0*f2;

        if f1*f1>4*f0*f2 then k0:=(-f1+sqrt(discr))/(2*f2) else k0:=-f1/(2*f2);
        if abs(k0)<PFAccuracy then k0:=(-f1-sqrt(discr))/(2*f2);
        k1:=e1+e2*k0;

        {Tabl:}
        tabl.w[j]:=w; tabl.k0[j]:=k0; tabl.k1[j]:=k1;
        if k0<>0 then tabl.k2[j]:=alfa*k1*k1/k0 else tabl.k2[j]:=0;
        w:=w+dw
      end;
    tabl.n:=nw
end;
{END - D-razbienie po n i m ------------------------------------------------------}

{BEGIN - Simou ------------------------------------------------------}
{procedure Simou(dt:real;n:integer;h:matr);}
procedure Simou(h : TablPK; inp : real; var PF : PFArray; var PFn : byte; var m,s : PFVector);
{h[0..n-1] - pereh krivaya, inp - vhodnoi signal; PF[1..9] - peredat. funkcii; m[0..4] - momenty, s[1..5] - ploschadi}
type matrice=array[0..5] of real;
var a,b : matrice;
    iter,k,i,mm,nn : integer;
    sum,s1,s2,fact,dt,kk : real;
begin
 dt:=h.t[1]-h.t[0]; h.n:=h.n-1; kk:=h.y[h.n];
 {Normirovanie:}
 if kk<>0 then for i:=0 to h.n do h.y[i]:=h.y[i]/kk else exit;
 {Raschet momentov (5 momentov):}
 for k:=0 to 4 do
 begin
  { Вычисление факторила }
  fact:=1;
  for i:=1 to k do  fact:=i*fact;
  { Вычисление моментов }
  s1:=0; s2:=0;
  if k=0 then
   for i:=1 to round(h.n/2) do
    begin
     s1:=s1+(1-h.y[2*i-1]);
     s2:=s2+(1-h.y[2*i]);
    end
  else
   for i:=1 to round(h.n/2) do
    begin
     s1:=s1+pow(k,-(2*i-1)*dt)*(1-h.y[2*i-1]);
     s2:=s2+pow(k,-(2*i)*dt)*(1-h.y[2*i]);
    end;
  if k=0 then m[k]:=(dt/fact/3)*((4*s1+2*s2+2-h.y[0]-h.y[h.n]))
  else m[k]:=(dt/fact/3)*(4*s1+2*s2+pow(k,-h.n*dt)*(1-h.y[h.n]));
  {writeln('m[',k,']=',m[k]:10:6);}
 end;
 { Вычисление площадей Симою }
 s[1]:=m[0];
 for k:=2 to 5 do
  begin
   sum:=0;
    for i:=0 to k-2 do sum:=sum+m[i]*s[k-1-i];
   s[k]:=m[k-1]+sum;
  end;
  {writeln('Вычисленные значения площадей Симою:');
  for k:=1 to 5 do
   writeln('s[',k,']=',s[k]:10:6);
  readln;}
  { Вычисление коофициентов передаточной функции }
  iter:=0;
  for mm:=0 to 2 do
   begin
    for i:=1 to 5 do  begin  a[i]:=0; b[i]:=0; end;
    a[0]:=1; b[0]:=1; s[0]:=1;
    for nn:=mm+1 to 5-mm do
     begin
      if mm=1 then b[mm]:=-s[nn+1]/s[nn];
      if mm=2 then
        begin
          b[mm]:=(sqr(s[nn+1])-s[nn]*s[nn+2])/(sqr(s[nn])-s[nn-1]*s[nn+1]);
          b[mm-1]:=(-s[nn+1]-b[mm]*s[nn-1])/s[nn];
        end;
    for k:=1 to nn do
       begin
        sum:=0;
        for i:=1 to k-1 do
          sum:=sum+b[i]*s[k-i];
        a[k]:=b[k]+s[k]+sum;
       end;
   iter:=iter+1; PFn:=iter;
   {clrscr;
   writeln('Вычисление коофициентов передаточной функции');
   writeln(iter,'-ый вариант при m=',mm,' n=',nn);}
   PF[iter].n:=nn; PF[iter].m:=mm; PF[iter].k:=kk/inp; PF[iter].tau:=0;
   for i:=0 to mm do PF[iter].b[i]:=b[i];
     {writeln('b[',i,']=',b[i]:7:3);
    writeln;}
    for i:=0 to nn do PF[iter].a[i]:=a[i];
     {writeln('a[',i,']=',a[i]:7:3);
   write('Для продолжения нажмите любую клавишу...');
    readln;}
   end;
   end;
   h.n:=h.n+1
end{Simou};
{END - Simou ------------------------------------------------------}





{BEGIN - Pokazateli kachestva ------------------------------------------------------}
procedure FindPPK(PK : TablPK; x : real; var PPK : TypePPK; var err : byte);
{Poisk pryamyh pokazateley. err=1 - process ne ustanovilsya, err=2 - shag po vremini <=0}
  var i : integer; y, yold, y1, y2, tr, y1max, y3max : real; poisk3 : boolean;
  begin
    err:=0;
    if PK.t[2]-PK.t[1]<PFAccuracy then begin err:=2; exit end;
    if abs(PK.y[PK.n]-PK.y[PK.n-1])<abs(0.05*PK.y[PK.n]) then
        PPK.yust:=PK.y[PK.n] else begin err:=1; exit end;
    PPK.est:=x-PPK.yust;
    y1max:=PPK.yust; y3max:=y1max; yold:=PK.y[0]; poisk3:=false;
    for i:=0 to PK.n do
      begin
      	y:=PK.y[i];
      	if (y>y1max) and (not poisk3) then y1max:=y; {Poisk A1}
      	if (y1max>PPK.yust) and (yold>y) then poisk3:=true; {Poisk A1 zavershen}
      	if (poisk3) and (yold<y) and (y>y3max) then y3max:=y; {Poisk A3}
      	yold:=y
      end;
    PPK.A1:=y1max-PPK.yust; PPK.A3:=y3max-PPK.yust;
    if (not poisk3) or (abs(PPK.A1)<PFAccuracy) then PPK.psi:=1 else PPK.psi:=1-PPK.A3/PPK.A1;
    if abs(PPK.yust)<PFAccuracy then PPK.sigma:=PPK.A1/PPK.yust else PPK.sigma:=PPK.A1;
    {Poisk tr:}
    if abs(PPK.yust)<PFAccuracy then {y -> 0; y1 - nijnyaya granitsa trubki, y2 - verhnyaya}
      begin y1:=-0.05*PPK.A1; y2:=-y1 end
      else begin y1:=0.95*PPK.yust; y2:=1.05*PPK.Yust end;
    for i:=0 to PK.n do
        if (PK.y[i]<=y1) or (PK.y[i]>=y2) then tr:=PK.t[i];
    PPK.tr:=tr
  end;

procedure FindKPK(pf : TransferFunction; var KPK : TypeKPK; var err : byte);
{Poisk kornevyh pokazateley. err=1 - net korney}
  var l, i : byte; x : PFComplexVector; m : real;
  begin
    KPK.n:=1000000; KPK.m:=1000000; KPK.ust:=true; err:=0;
    {Korni:}
    Roots(pf.a,pf.n,x,l);
    if l=0 then begin err:=1; exit end;
    {Analiz korney:}
    for i:=0 to l do
      begin
      	if -x[i].re<KPK.n then KPK.n:=-x[i].re;
      	if abs(x[i].im)>PFAccuracy then m:=abs(x[i].re/x[i].im) else m:=1000000;
      	if KPK.m>m then KPK.m:=m
      end;
    if KPK.n<0 then KPK.ust:=false
  end;

procedure FindFPK(tablWraz, tablFz, tablFe : TablAFH; var FPK : TypeFPK; var err : byte);
{Poisk chastotnyh pokazateley. tablWraz, tablFz, tablFe - AFH Wraz, Fz i Fe}
{err=1 - Az(w) rastet, no ne dostigaet M; err=2 - Ae(w) rastet, no ne dostigaet ME; err=3 - shag po chastote <=0}
  var i, j : integer; m, w : real;
  begin
    err:=0;
    if tablWraz.w[2]-tablWraz.w[1]<PFAccuracy then begin err:=3; exit end;
    m:=tablFz.A[0]; w:=tablFz.w[0];
    for i:=0 to tablFz.n do if m<tablFz.A[i] then
      begin
      	m:=tablFz.A[i]; w:=tablFz.w[i];
      	if i=tablFz.n then err:=1
      end;
    FPK.M:=m; FPK.wr:=w;
    m:=tablFe.A[0]; w:=tablFe.w[0];
    for i:=0 to tablFe.n do if m<tablFe.A[i] then 
      begin 
      	m:=tablFe.A[i]; w:=tablFe.w[i];
      	if i=tablFe.n then err:=2
      end;
    FPK.ME:=m; FPK.we:=w;
    m:=tablWraz.A[0]; j:=0;
    for i:=tablWraz.n downto 1 do 
      begin
      	if (tablWraz.re[i]<0) and (tablWraz.im[i]*tablWraz.im[i-1]<0) then m:=-tablWraz.re[i];
      	if (tablWraz.a[i]<1) and (tablWraz.a[i-1]>1) then j:=i
      end;
    FPK.dA:=m; FPK.dfi:=tablWraz.fi[j]-180
  end;
{END - Pokazateli kachestva ------------------------------------------------------}




{BEGIN - Perehodnye processy v ASR ------------------------------------------------------}
procedure ASRProcesses(wr,wob : TransferFunction; y0 : PFVector; y0n : byte; x : real; t0, dt : real; nt : byte;
var yt, et, ut : TablPK);
{Raschet perehodnyh processov v ASR: wr, wob - PF regulyatora i OU, y0 - vector nach usloviy dlya OU}
  var r : Block2; ob : BlockN; tau : BlockT; i : byte; y, e, yy, u, t : real;
  begin
    if (nt=0) or (dt<=0) then begin et.n:=0; ut.n:=0; yt.n:=0; exit end;
    TimeDiscr:=dt;
    r.init(wr.k,wr.a[0],wr.a[1],wr.a[2],wr.b[0],wr.b[1],wr.b[2],0,0,0,0);
    ob.zeros; ob.k:=wob.k;
    for i:=0 to wob.m do ob.b[i]:=wob.b[i];
    for i:=0 to wob.n do ob.a[i]:=wob.a[i];
    for i:=0 to maxn do begin ob.xt[i]:=0; ob.yt[i]:=y0[0] end;
    ob.init;
    yt.n:=nt; et.n:=nt; ut.n:=nt;
    tau.init(wob.tau,y0[0]);
    y:=y0[0]; u:=0;  t:=t0;
    for i:=0 to nt do
      begin
      	e:=x-y;
      	u:=r.run(e);
        yt.y[i]:=y; et.y[i]:=e; ut.y[i]:=u;
      	yy:=ob.run(u);
      	y:=tau.run(yy);
        yt.t[i]:=t; et.t[i]:=t; ut.t[i]:=t; t:=t+dt
      end
  end;
{END - Perehodnye processy v ASR ------------------------------------------------------}





{BEGIN - D-razbienie po M i ME ------------------------------------------------------}
procedure dmme(ws : TransferFunction; m : real; typ : byte; alfa : real; w0, dw : real; nw : byte;
var tabl : TablD);
{Krivaya D-razb dlya:}
{typ=0 - raschet dlya M, typ=1 - dlya ME}
  var
    i, j, k, kwn, kmn : byte;
    mm, wm : dmatr;
    w, c, kw, km, y0w, y0m, dk0, dk1, k00, k10, k0itog, k1itog, wmain : real;
    k0w, k1w, k0m, k1m : array[1..4] of real;
    wdone : boolean;
    wr : TransferFunction;
  label Loop;

  procedure f(k0,k1 : real; var m,w : real);
    {Procedura rascheta znacheniy po zadannym parametram}
    var tablWraz, tablD, tablFz, tablFe, tablFv : TablAFH; FPK : TypeFPK; err : byte; dww : real;
    begin
      wr.b[0]:=k0; wr.b[1]:=k1;
      if abs(k0)>PFAccuracy then wr.b[2]:=alfa*k1*k1/k0 else wr.b[2]:=1000000;
      case typ of
      0: begin {po M}
           dww:=0.2*wmain;
           repeat
             afhasr(wr,ws,0,0,dww,20,0,tablWraz,tablD,tablFz,tablFe,tablFv);
             FindFPK(tablWraz,tablFz,tablFe,FPK,err);
             dww:=dww*5
           until (err=0);
           m:=FPK.M; w:=FPK.wr
      	 end;
      1: begin {po ME}
           dww:=0.2*wmain;
           repeat
             afhasr(wr,ws,0,0,dww,20,0,tablWraz,tablD,tablFz,tablFe,tablFv);
             FindFPK(tablWraz,tablFz,tablFe,FPK,err);
             dww:=dww*5
           until (err=0);
           m:=FPK.ME; w:=FPK.we
      	 end
      end
    end;

  begin
    if ws.k<PFAccuracy then exit;
    if nw=0 then exit;
    dw:=abs(dw);
    dk0:=0.1/ws.k; dk1:=dk0; k00:=-0.4*maxmatr*dk0; k10:=k00;
    wr.k:=1; wr.tau:=0; wr.a[0]:=0; wr.a[1]:=1; wr.n:=1; wr.m:=2;
    wmain:=nw*dw;
    {Zapolnenie matricy:}
    for i:=0 to maxmatr do for j:=0 to maxmatr do f(k00+dk0*i,k10+dk1*j,mm[i,j],wm[i,j]);
    k:=1;
{    for k:=1 to nw do}
Loop: begin
      	w:=w0+(k-1)*dw;
      	wdone:=false;
      	for i:=0 to maxmatr-1 do
      	    for j:=0 to maxmatr-1 do
      	      begin
{MessageDlg(inttostr(i)+' '+inttostr(j)+' '+floattostr(w), mtInformation,[mbOk], 0);}
      	      	if ((wm[i,j]<w) and (wm[i+1,j]<w) and (wm[i,j+1]<w) and (wm[i+1,j+1]<w)) or
      	      	   ((wm[i,j]>w) and (wm[i+1,j]>w) and (wm[i,j+1]>w) and (wm[i+1,j+1]>w)) or
      	      	   ((mm[i,j]<m) and (mm[i+1,j]<m) and (mm[i,j+1]<m) and (mm[i+1,j+1]<m)) or
      	      	   ((mm[i,j]>m) and (mm[i+1,j]>m) and (mm[i,j+1]>m) and (mm[i+1,j+1]>m)) or
      	      	   (wdone) then
      	      	   begin {kletka ne podhodit}
      	      	   end
      	      	   else begin
      	      	     {Prosmotr chastot:}
      	      	     kwn:=0;
      	      	     {1-ya gran'}
      	      	     if wm[i+1,j]<>wm[i+1,j+1] then
      	      	       begin
      	      	       	 c:=(wm[i+1,j+1]-w)/(wm[i+1,j+1]-wm[i+1,j]);
      	      	       	 if (c>=0) and (c<=1) then
      	      	       	   begin
      	      	       	     kwn:=kwn+1; k0w[kwn]:=k00+dk0*(i+1); k1w[kwn]:=k10+dk1*(j+1)-c*dk1
      	      	       	   end
      	      	       end;
      	      	     {2-ya gran'}
      	      	     if wm[i,j]<>wm[i+1,j] then
      	      	       begin
      	      	       	 c:=(wm[i+1,j]-w)/(wm[i+1,j]-wm[i,j]);
      	      	       	 if (c>=0) and (c<=1) then
      	      	       	   begin
      	      	       	     kwn:=kwn+1; k1w[kwn]:=k10+dk1*j; k0w[kwn]:=k00+dk0*(i+1)-c*dk0
      	      	       	   end
      	      	       end;
      	      	     {3-ya gran'}
      	      	     if wm[i,j]<>wm[i,j+1] then
      	      	       begin
      	      	       	 c:=(wm[i,j+1]-w)/(wm[i,j+1]-wm[i,j]);
      	      	       	 if (c>=0) and (c<=1) then
      	      	       	   begin
      	      	       	     kwn:=kwn+1; k0w[kwn]:=k00+dk0*i; k1w[kwn]:=k10+dk1*(j+1)-c*dk1
      	      	       	   end
      	      	       end;
      	      	     {4-ya gran'}
      	      	     if wm[i,j+1]<>wm[i+1,j+1] then
      	      	       begin
      	      	       	 c:=(wm[i+1,j+1]-w)/(wm[i+1,j+1]-wm[i,j+1]);
      	      	       	 if (c>=0) and (c<=1) then
      	      	       	   begin
      	      	       	     kwn:=kwn+1; k1w[kwn]:=k10+dk1*(j+1); k0w[kwn]:=k00+dk0*(i+1)-c*dk0
      	      	       	   end
      	      	       end;
      	      	     {Raschet kw i y0w:}
      	      	     if kwn>=2 then
      	      	       begin
      	      	         if k1w[2]<>k1w[1] then kw:=(k0w[2]-k0w[1])/(k1w[2]-k1w[1]) else kw:=1000000;
                         y0w:=k0w[1]-kw*k1w[1]
      	      	       end;

      	      	     {Prosmotr M:}
      	      	     kmn:=0;
      	      	     {1-ya gran'}
      	      	     if mm[i+1,j]<>mm[i+1,j+1] then 
      	      	       begin
      	      	       	 c:=(mm[i+1,j+1]-m)/(mm[i+1,j+1]-mm[i+1,j]);
      	      	       	 if (c>=0) and (c<=1) then
      	      	       	   begin
      	      	       	     kmn:=kmn+1; k0m[kmn]:=k00+dk0*(i+1); k1m[kmn]:=k10+dk1*(j+1)-c*dk1
      	      	       	   end
      	      	       end;
      	      	     {2-ya gran'}
      	      	     if mm[i,j]<>mm[i+1,j] then
      	      	       begin
      	      	       	 c:=(mm[i+1,j]-m)/(mm[i+1,j]-mm[i,j]);
      	      	       	 if (c>=0) and (c<=1) then
      	      	       	   begin
      	      	       	     kmn:=kmn+1; k1m[kmn]:=k10+dk1*j; k0m[kmn]:=k00+dk0*(i+1)-c*dk0
      	      	       	   end
      	      	       end;
      	      	     {3-ya gran'}
      	      	     if mm[i,j]<>mm[i,j+1] then 
      	      	       begin
      	      	       	 c:=(mm[i,j+1]-m)/(mm[i,j+1]-mm[i,j]);
      	      	       	 if (c>=0) and (c<=1) then
      	      	       	   begin
      	      	       	     kmn:=kmn+1; k0m[kmn]:=k00+dk0*i; k1m[kmn]:=k10+dk1*(j+1)-c*dk1
      	      	       	   end
      	      	       end;
      	      	     {4-ya gran'}
      	      	     if mm[i,j+1]<>mm[i+1,j+1] then
      	      	       begin
      	      	       	 c:=(mm[i+1,j+1]-m)/(mm[i+1,j+1]-mm[i,j+1]);
      	      	       	 if (c>=0) and (c<=1) then
      	      	       	   begin
      	      	       	     kmn:=kmn+1; k1m[kmn]:=k10+dk1*(j+1); k0m[kmn]:=k00+dk0*(i+1)-c*dk0
      	      	       	   end
      	      	       end;
      	      	     {Raschet km i y0m:}
      	      	     if kmn>=2 then
      	      	       begin
      	      	         km:=(k0m[2]-k0m[1])/(k1m[2]-k1m[1]);
                         y0m:=k0m[1]-km*k1m[1]
      	      	       end;

      	      	     {Koordinata tochki peresecheniya:}
      	      	     if (kwn>=2) and (kmn>=2) then
      	      	       begin
      	      	       	 if kw<>km then k1itog:=(y0m-y0w)/(kw-km) else k1itog:=1000000;
      	      	       	 k0itog:=kw*k1itog+y0w;
      	      	       	 if (k0itog>k00+dk0*i) and (k0itog<k00+dk0*(i+1)) and
      	      	       	    (k1itog>k10+dk1*j) and (k1itog<k10+dk1*(j+1)) then
      	      	       	    begin 
      	      	       	      tabl.w[k]:=w;
{MessageDlg('hdfjh jfhjhjhj hjhjhjfdhjhjhjhjhjdfhgj hfds '+floattostr(w), mtInformation,[mbOk], 0);}
      	      	       	      tabl.k1[k]:=k1itog; tabl.k0[k]:=k0itog;
      	      	       	      if abs(k0itog)>PFAccuracy then tabl.k2[k]:=alfa*k1itog*k1itog/k0itog else tabl.k2[k]:=1000000;
      	      	       	      wdone:=true
      	      	       	    end
      	      	       end
      	      	   end {else}
      	      end; {i,j}
      	if not wdone then {Pererazbienie ploskosti}
      	  begin
            dk0:=2*dk0; dk1:=2*dk1; k00:=2*k00; k10:=2*k10;
{MessageDlg(floattostr(k00)+' '+floattostr(k10)+' '+floattostr(dk0)+' '+floattostr(dk1), mtInformation,[mbOk], 0);}
            {Zapolnenie matricy:}
            for i:=0 to maxmatr do for j:=0 to maxmatr do f(k00+dk0*i,k10+dk1*j,mm[i,j],wm[i,j]);
            k:=k-1
      	  end
      end; {k}
    k:=k+1;
    if k<=nw then goto Loop;
    tabl.n:=nw
  end;
{END - D-razbienie po M i ME ------------------------------------------------------}

end.





{FOR DEBUGGING ------------------------------------------------------}
(*
var m,n,l,j,PFn:byte; a,b,c:PFVector; x,km:PFComplexVector; xx:complex; i:byte;
  compl, krat, nom : PFByteVector; h:TablPK; PF:PFArray;
begin
  writeln;
  GetParmsFromStr('-12.3s5 256.78-95 -12..5 adassa 23.6',a,n);
  writeln(n);
  for i:=0 to n-1 do writeln(a[i]:8:4);
end.      *)
(*
  for i:=0 to maxn do h[i].t:=i;
  h[0].y:=0; h[1].y:=2.2; h[2].y:=6; h[3].y:=9.2; h[4].y:=11.6;
  h[5].y:=13.8; h[6].y:=15.7; h[7].y:=17.5; h[8].y:=19.1; h[9].y:=20.4;
  h[10].y:=21.3; h[11].y:=21.9; h[12].y:=22.3; h[13].y:=22.5; h[14].y:=22.6;
  n:=14;
  Simou(h,n,PF,PFn,a,b);
  writeln('n=',n);
  for i:=0 to n do write(' ',i,' ',h[i].y:8:4); writeln;
  for i:=1 to 5 do writeln(a[i-1]:8:4,b[i]:8:4); writeln;
  writeln('PFn=',PFn);
  readln;
  for i:=1 to PFn do
    begin
      writeln('k,tau=',PF[i].k:8:4,PF[i].tau:8:4);
      for j:=0 to PF[i].m do write(' ',PF[i].b[j]:8:4); writeln;
      for j:=0 to PF[i].n do write(' ',PF[i].a[j]:8:4); writeln;
    end;

end.
  for i:=0 to maxn do a[i]:=0;
{  a[4]:=1; a[3]:=-2; a[2]:=-17; a[1]:=18; a[0]:=72; m:=4;}
{  a[4]:=1; a[3]:=1; a[2]:=1; a[1]:=1; a[0]:=1; m:=2;}
  a[4]:=2.5; a[3]:=5.75; a[2]:=4.25; a[1]:=5.35; a[0]:=0; m:=4;
{  a[3]:=2.5; a[2]:=5.75; a[1]:=4.25; a[0]:=5.35; m:=3;}
{  a[2]:=1; a[1]:=5; a[0]:=6; m:=2;}
{  a[2]:=1; a[1]:=2; a[0]:=1; m:=2;}
{   a[5]:=1; a[4]:=5; a[3]:=6; a[2]:=0; a[1]:=0; a[0]:=0; m:=5;}
{   a[4]:=1; a[3]:=-2; a[2]:=-3; a[1]:=4; a[0]:=4; m:=4;}
  Roots(a,m,x,n);
  writeln('n = ',n:3);
  for j:=1 to n do writeln('x[',j,'] = ',x[j].stroka);

  b[0]:=4.35;
{  LaplasM(1,a,m,b,0,0,x,km,l,compl,krat,nom);}
  writeln('l=',l);
  for i:=1 to l do
      writeln('x,m,c,k,n=',x[i].stroka,km[i].stroka,compl[i],krat[i],nom[i]);
end.*)