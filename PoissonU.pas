{$N+,R+,J+}
{ define STANDALONE}
{$ifdef STANDALONE}
{ define USE_ARRAY}
program Poisson_; {generates random numbers in Poisson distribution}
uses OPCRT,GRAPH,DOS,MRNG;
{$else}
unit PoissonU;

INTERFACE
   uses MRNG,SYSUTILS;
{$ifdef NDEBUG}
type double=real;
{$endif NDEBUG}
function Factorial(n :integer):double;
function Poisson(x :double; n :integer):double;
function OldFastPoissonDist(xm :double):integer; {returns a random number in Poisson distr w/ mean xm}
function FastPoissonDist(xm :double):integer; {returns a random number in Poisson distr w/ mean xm}
function PoissonDist(xm :double; max :integer):integer; {returns a random number in Poisson distr w/ mean xm}
function IntegerPoissonDist(xm :integer; max :integer):integer; {returns a random number in Poisson distr w/ mean xm}
function ApproxPoissonDist(xm :double; max :integer):integer; {returns a random number in Poisson distr w/ mean xm}
function OldPoissonFromSortedArray(xm :double):integer; {returns a random number in Poisson distr w/ mean xm}
function PoissonFromSortedArray(xm :double):integer; {returns a random number in Poisson distr w/ mean xm}
(* function TableAssistedPoissonDist(xm :double):integer; {returns a random number in Poisson distr w/ mean xm}
   I don't trust this function any more - produces biased results. *)

IMPLEMENTATION
{$endif}

{type double=real;}
const maxarg=40;
      maxres=60;
      maxintarg=70;
      maxintres=107;
      maxsqr=127;
      lastx:double=-1;
      maxarray=10000000;
      maxfact=20;

type realarray=array[0..maxres] of double;
     largearray=array[0..maxsqr] of realarray;
     xrealarray=array[0..maxintres] of double;
     pfintarray=array[1..maxintarg] of xrealarray;
     bigarray=array[1..maxarray] of double;

var pf  :^pfintarray;
    lpa :^largearray;
    rnd :^bigarray;
    p   :array[0..maxintres] of double;
    fact         :array[0..maxfact] of double;
    fs,top,ii    :longint;
    invfs        :double;
    randfile     :file of double;

{$ifdef STANDALONE}
procedure VGAMode;
          var driver, mode: integer;
          begin
          driver:=VGA;
          mode:=VGAHI;
          InitGraph(driver,mode,'C:\BP\BGI');
          if (grerror<>0) then writeln(grerror,': ',GraphErrorMsg(grerror));
          SetTextStyle(smallfont,horizdir,4);
          SetFillStyle(solidfill,cyan);
          end;
{$endif STANDALONE}

function Factorial(n :integer):double;
         var i      :integer;
         begin
         if (n<=maxfact) then Factorial:=fact[n];
         {if (n>170) then begin
            writeln('Overflow: Factorial called for argument >170');
            halt;
            end;}
         result:=1;
         for i:=1 to n do result:=result*i;
         end;

function Poisson(x :double; n :integer):double;
         begin
         {writeln(n:4,x:10);}
         if (n>170) then Poisson:=exp(n*ln(x-n/exp(1))-x)/sqrt(2*pi*n)
         else Poisson:=exp(n*ln(x)-x)/factorial(n);
         end;

procedure InitializeLargeArray;
          var i,j :integer;
              x   :double;
          begin
          fact[0]:=1;
          for i:=1 to maxfact do fact[i]:=fact[pred(i)]*i;
          new(lpa);
          lpa^[0,0]:=1;
          for j:=1 to 60 do begin
            lpa^[0,j]:=0;
            end;
          for i:=1 to maxsqr do begin
            x:=sqr(i)/512;
            for j:=0 to 60 do begin
              lpa^[i,j]:=Poisson(x,j);
              end;
            end;
          end;

procedure ComputePF;
          var i,j :integer;
          begin
          new(pf);
          for i:=1 to maxintarg do
             for j:=0 to maxintres do
                pf^[i,j]:=Poisson(i,j);
          end;

function GammLn(xx :double):double; {returns ln(Gamma(xx))}
         const stp=2.50662827465;   {Numerical Recipes, p 704}
               half=0.5; one=1.0; fpf=5.5;
               cof :array[1..6] of double=(76.18009173,-86.50532033,24.0109822,
                                           -1.231739516,0.12058003E-2, -0.536382E-5);
         var x,tmp,ser :double; j:integer;
         begin
         x:=xx-one; tmp:=x+fpf; tmp:=(x+half)*ln(tmp)-tmp; ser:=one;
         for j:=1 to 6 do begin
           x:=x+one; ser:=ser+cof[j]/x; end;
         GammLn:=tmp+ln(stp*ser);
         end;


procedure CreatePoissonArray(xm :double; max :integer);
          var i :integer;
          begin
          for i:=0 to max do p[i]:=Poisson(xm,i);
          lastx:=xm;
          end;

function OldFastPoissonDist(xm :double):integer; {returns a random number in Poisson distr w/ mean xm}
         var r2xm              :integer;
             cumprob,newprob   :double;
         begin
         result:=round(xm);
         newprob:=Poisson(xm,result);
         if (mrandom<newprob) then exit;
         r2xm:=2*result; cumprob:=1;
         repeat
           cumprob:=cumprob-newprob;
           result:=pred(r2xm-result);
           if (result>=0) then begin
             newprob:=Poisson(xm,result);
             if (mrandom<newprob/cumprob) then exit;
             cumprob:=cumprob-newprob;
             end;
           result:=r2xm-result;
           newprob:=Poisson(xm,result);
           if {(result>169) or} (mrandom<newprob/cumprob) then exit;
         until false;
         end;

function FastPoissonDist(xm :double):integer; {returns a random number in Poisson distr w/ mean xm}
         var lores,hires            :integer;
             cumprob,newprob        :double;
         begin
         lores:=round(xm);
         newprob:=Poisson(xm,lores);
         if (mrandom<newprob) then begin FastPoissonDist:=lores; exit; end;
         cumprob:=1; hires:=lores;
         while (lores>0) do begin
           dec(lores);
           cumprob:=cumprob-newprob;
           newprob:=Poisson(xm,lores);
           if (mrandom<newprob/cumprob) then begin FastPoissonDist:=lores; exit; end;
           cumprob:=cumprob-newprob;
           inc(hires);
           newprob:=Poisson(xm,hires);
           if (mrandom<newprob/cumprob) then begin FastPoissonDist:=hires; exit; end;
           end;
         repeat
           cumprob:=cumprob-newprob;
           inc(hires);
           newprob:=Poisson(xm,hires);
           if (mrandom<newprob/cumprob) then begin FastPoissonDist:=hires; exit; end;
         until false;
         end;

function ApproxPoissonDist(xm :double; max :integer):integer; {returns a random number in Poisson distr w/ mean xm}
         var i :integer;
             p :double;
         begin
         p:=xm/max;
         result:=0;
         for i:=1 to max do if (mrandom<p) then inc(result);
         end;

function SlowPoissonDist(xm :double; max :integer):integer; {returns a random number in Poisson distr w/ mean xm}
         begin
         {if (xm<>lastx) then} CreatePoissonArray(xm,max);
         repeat result:=mrandint(0,pred(max)) until (mrandom<p[result]);
         end;

function PoissonDist(xm :double; max :integer):integer; {returns a random number in Poisson distr w/ mean xm}
          {Divides the domain in half as many times as necessary}
         var i,vhi,vlo,vmid :integer;
             s1,s2       :double;
         begin
         {if (xm<>lastx) then} CreatePoissonArray(xm,max);
         vhi:=max; vlo:=0;
         vmid:=round(xm)+2;
         repeat
           s1:=0; s2:=0;
           for i:=vlo to vmid do s1:=s1+p[i];
           for i:=succ(vmid) to vhi do s2:=s2+p[i];
           if (mrandom<s1/(s1+s2)) then {low half of range chosen}
             vhi:=vmid
           else
             vlo:=succ(vmid);
           vmid:=(vhi+vlo) div 2;
         until (vlo=vhi);
         PoissonDist:=vlo;
         end;

function TableAssistedPoissonDist(xm :double):integer; {returns a random number in Poisson distr w/ mean xm}
          {Divides the domain in half as many times as necessary}
         var i,vhi,vlo,vmid,max :integer;
             s1,s2       :double;
             p           :^realarray;
         begin
         i:=round(sqrt(512*xm));
         max:=round(xm)+i div 10;
         if (max>60) then max:=60;
         p:=@lpa^[i];
         vhi:=max; vlo:=0;
         vmid:=round(xm)+2;
         repeat
           s1:=0; s2:=0;
           for i:=vlo to vmid do s1:=s1+p^[i];
           for i:=succ(vmid) to vhi do s2:=s2+p^[i];
           if (mrandom<s1/(s1+s2)) then {low half of range chosen}
             vhi:=vmid
           else
             vlo:=succ(vmid);
           vmid:=(vhi+vlo) div 2;
         until (vlo=vhi);
         TableAssistedPoissonDist:=vlo;
         end;

function IntegerPoissonDist(xm :integer; max :integer):integer; {returns a random number in Poisson distr w/ mean xm}
          {Quick shortcut version of PoissonDist above, xm an integer}
         var i,vhi,vlo,vmid :integer;
             s1,s2       :double;
             p           :^xrealarray;
         begin
         if (xm>maxintarg) then begin
           IntegerPoissonDist:=ApproxPoissonDist(xm,max);
           exit;
           end;
         p:=@pf^[xm];
         vhi:=max; vlo:=0;
         vmid:=round(xm)+2;
         repeat
           s1:=0; s2:=0;
           for i:=vlo to vmid do s1:=s1+p^[i];
           for i:=succ(vmid) to vhi do s2:=s2+p^[i];
{if (s1+s2=0) then
  writeln('ERROR IN INTEGER POISSON: ',xm:4,max:4);}
           if (mrandom<s1/(s1+s2)) then {low half of range chosen}
             vhi:=vmid
           else
             vlo:=succ(vmid);
           vmid:=(vhi+vlo) div 2;
         until (vlo=vhi);
         IntegerPoissonDist:=vlo;
         end;
(*
function PoissonDev(xm :double;max :integer):double; {returns a random number in Poisson distr w/ mean xm}
         const gloldm :double=-1;   {Numerical Recipes, p 717}
               one :double=1;
         var glsq,glalxm,glg  :double;
             em,t,y           :double;
         begin
         if (xm<12) then begin
           if (xm<gloldm) then begin gloldm:=xm; glg:=exp(-xm); end;
           em:=-1;
           t:=1;
           repeat
             em:=em+one; t:=t*mrandint(0,pred(max));
           until (t<=glg);
           end
         else begin
           if (xm<>gloldm) then begin
             gloldm:=xm; glsq:=sqrt(2*xm); glalxm:=ln(xm);
             glg:=xm*glalxm-gammln(xm+one);
             end;
           repeat
             repeat
               y:=y*mrandint(0,pred(max)); y:=sin(y)/cos(y); em:=glsq*y*xm;
             until (em>=0);
             em:=trunc(em);
             t:=0.9*(one+sqr(y))*exp(em*glalxm-gammln(em+one)-glg);
           until (mrandint(0,pred(max))<=t);
           end;
         PoissonDev:=em;
         end;
*)
function FindPos(r :double):longint;
         var jump,rslt :longint;
             q         :double;
         begin
         rslt:=round(fs*r);
         seek(randfile,rslt); read(randfile,q);
         jump:=5+round(abs(q-r)*5*fs);
         repeat
           if (q>r) then dec(rslt,jump) else inc(rslt,jump);
           jump:=succ(jump) div 2;
           seek(randfile,rslt); read(randfile,q);
         until (jump=1);
         if (q>r) then
           repeat dec(rslt); seek(randfile,rslt); read(randfile,q); until (q<r)
         else begin repeat inc(rslt); read(randfile,q); until (q>r); dec(rslt); end;
         FindPos:=rslt;
         end;

function PoissonFromSortedFile(xm :double):integer; {returns a random number in Poisson distr w/ mean xm}
          var r1,r2 :double;
          begin
          r1:=mrandom;
          r2:=r1+invfs*xm;
          if (r2>1) then begin r2:=r1; r1:=r2-invfs*xm; end;
          PoissonFromSortedFile:=FindPos(r2)-FindPos(r1);
          end;

function FasterPoissonFromSortedFile(xm :double):integer; {returns a random number in Poisson distr w/ mean xm}
          var r,r1,r2 :double;
              p,p1      :longint;
          begin
          p1:=mrandint(0,pred(top));
          seek(randfile,p1); read(randfile,r1);
          r2:=r1+invfs*xm;
          p:=-1;
          repeat read(randfile,r); inc(p) until (r>r2);
          FasterPoissonFromSortedFile:=p;
          end;

function PoissonFromSortedArray(xm :double):integer; {returns a random number in Poisson distr w/ mean xm}
          var r1,r2 :double;
              p,p1      :longint;
          begin
          repeat
            r1:=mrandom
          until (r1<0.9995);
          r2:=r1+invfs*xm;
          p:=round(fs*r1); {Estimate where r1 is in randfile}
          p:=p+round(fs*(r1-rnd^[p])); {refine search for r1, twice}
          p:=p+round(fs*(r1-rnd^[p]));
          if (rnd^[p]<r1) then begin
            repeat inc(p) until (rnd^[p]>r1); p1:=p;
            while (rnd^[p]<r2) do inc(p);
            PoissonFromSortedArray:=p-p1;
            end
          else if (rnd^[p]<r2) then begin
            p1:=p;
            repeat dec(p1) until (rnd^[p1]<r1);
            repeat inc(p) until (rnd^[p]>r2); dec(p);
            PoissonFromSortedArray:=p-p1;
            end
          else begin
            repeat dec(p) until (rnd^[p]<r2); p1:=p;
            while (rnd^[p1]>r1) do dec(p1);
            PoissonFromSortedArray:=p-p1;
            end;
          end;

function OldPoissonFromSortedArray(xm :double):integer; {returns a random number in Poisson distr w/ mean xm}
          var r1,r2 :double;
              p,p1      :longint;
          begin
          {p1:=mrandint(0,pred(top));}
          p1:=random(top);
          r1:=rnd^[p1];
          r2:=r1+invfs*xm;
          p:=p1;
          repeat inc(p) until (rnd^[p]>r2); dec(p);
          OldPoissonFromSortedArray:=p-p1;
          end;

{$ifdef STANDALONE}
var i,x,y,n :integer;
    h1,m1,s1,sh1,h2,m2,s2,sh2 :word;
    r,sum   :double;
    s       :string;
    hist    :array[0..255] of integer;
    cumtime :array[0..10] of integer;
    name    :array[0..10] of string;

begin
assign(randfile,{'C:\AGNERFOG\DATA\}'E:RANDNUM.DAT'); reset(randfile);
fs:=filesize(randfile); invfs:=1/fs; top:=round(0.99*fs); {jumpstart:=5*round(sqrt(fs));}
new(rnd); seek(randfile,random(top)); for i:=1 to maxarray do read(randfile,rnd^[i]);

i:=30;
{repeat inc(i); writeln(i,factorial(i)); until false;}
randomize;
InitializeLargeArray;
ComputePF;
{for i:=1 to 500 do write(PoissonDist(6,15):8:3);
writeln(factorial(20):28:2,factorial(19):28:2,factorial(20)/factorial(19):8:2);}
VGAMode;
sum:=0;
for n:=0 to 50 do begin
   hist[n]:=0;
   str(n:2,s);
   x:=12*n;
   y:=479-4*round(500*Poisson(5,n));
   sum:=sum+Poisson(10,n);
   Bar(x,y,x+12,479);
   if (n mod 2 =0) then OutTextXY(x,465,s)
   {Line(x,y,x,479);}
   end;

SetColor(yellow);
GetTime(h1,m1,s1,sh1);
for i:=1 to 30000 do inc(hist[PoissonDist(5,15)]);
GetTime(h2,m2,s2,sh2); cumtime[0]:=sh2-sh1 + 100*(s2-s1) + 6000*(m2-m1); name[0]:='PoissonDist';
for n:=0 to 50 do begin line(7+12*n,479,7+12*n,479-(hist[n] div 15)); hist[n]:=0; end;

SetColor(blue);
GetTime(h1,m1,s1,sh1);
for i:=1 to 30000 do inc(hist[IntegerPoissonDist(5,15)]);
GetTime(h2,m2,s2,sh2); cumtime[1]:=sh2-sh1 + 100*(s2-s1) + 6000*(m2-m1); name[1]:='IntegerPoissonDist';
for n:=0 to 50 do begin line(9+12*n,479,9+12*n,479-(hist[n] div 15)); hist[n]:=0; end;

SetColor(lightred);
GetTime(h1,m1,s1,sh1);
for i:=1 to 30000 do inc(hist[SlowPoissonDist(5,15)]);
GetTime(h2,m2,s2,sh2); cumtime[2]:=sh2-sh1 + 100*(s2-s1) + 6000*(m2-m1); name[2]:='SlowPoissonDist';
for n:=0 to 50 do begin line(5+12*n,479,5+12*n,479-(hist[n] div 15)); hist[n]:=0; end;

SetColor(darkgray);
GetTime(h1,m1,s1,sh1);
{for i:=1 to 30000 do inc(hist[TableAssistedPoissonDist(5)]);}
{for i:=1 to 30000 do inc(hist[FasterPoissonFromSortedFile(5)]);}
for i:=1 to 30000 do inc(hist[PoissonFromSortedArray(5)]);
GetTime(h2,m2,s2,sh2); cumtime[5]:=sh2-sh1 + 100*(s2-s1) + 6000*(m2-m1); name[5]:='PoissonDistFromSorted';
for n:=0 to 50 do begin line(13+12*n,479,13+12*n,479-(hist[n] div 15)); hist[n]:=0; end;SetColor(lightmagenta);

SetColor(white);
GetTime(h1,m1,s1,sh1);
for i:=1 to 30000 do inc(hist[FastPoissonDist(5)]);
GetTime(h2,m2,s2,sh2); cumtime[3]:=sh2-sh1 + 100*(s2-s1) + 6000*(m2-m1); name[3]:='FastPoissonDist';
for n:=0 to 50 do begin line(11+12*n,479,11+12*n,479-(hist[n] div 15)); hist[n]:=0; end;

SetColor(lightgreen);
GetTime(h1,m1,s1,sh1);
for i:=1 to 30000 do inc(hist[ApproxPoissonDist(5,45)]);
GetTime(h2,m2,s2,sh2); cumtime[4]:=sh2-sh1 + 100*(s2-s1) + 6000*(m2-m1); name[4]:='ApproxPoissonDist';
for n:=0 to 50 do begin line(3+12*n,479,3+12*n,479-(hist[n] div 15)); hist[n]:=0; end;


repeat until keypressed;
CloseGraph;
writeln(sum);
for i:=0 to 5 do writeln(cumtime[i]:5,'  ',name[i]);
{$else STANDALONE}
begin
MRandSeed(round(maxint*frac(time)));
InitializeLargeArray;
ComputePF;
{$ifdef USE_ARRAY}
assign(randfile,'C:\RANDOM\MRANDNUM.100'); reset(randfile);
fs:=filesize(randfile); invfs:=1/fs; top:=round(0.99*fs);
new(rnd);
if (maxarray>fs) then writeln('Array ',maxarray,' too large.');
if (maxarray<fs) then seek(randfile,random(top));
for ii:=1 to maxarray do read(randfile,rnd^[ii]);
{$endif USE_ARRAY}
{$endif STANDALONE}
end.
