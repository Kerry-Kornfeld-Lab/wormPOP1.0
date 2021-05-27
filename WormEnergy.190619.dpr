program WormEnergy;
   {Kornfeld's C elegans experiment in population dynamics}
{$apptype CONSOLE}
{$N+,B-,O-,J+,R+,E-,X+,V+}
{$R+}

{ define NDEBUG}
{ define RDEBUG}
{$define INDOUTPUT}
{ define PHEROMONE}
{ define VARIABILITY}  {This option makes appetite constant and min lean mass into individual variables, for each worm, distributed about an input mean}

uses
  SysUtils,
  DOSSTUFF,
{$ifdef VARIABILITY}
  GAUSSDIS,
{$endif VARIABILITY}
  POISSONU;

const vxstr='WormEnergy vx190619 ';

const crlf=#13#10; maxlife=499; {periods}  one:double=1.0; half:double=0.5; quarter:double=0.25;  eighth=0.125; tenth:double=0.1; tiny:double=1.0E-12; small=0.01; ln2=0.69314718056; infinity=1.0E30; sq2=1.41421356;  twothirds=2/3; threehalves=1.5;
      {$ifdef INDOUTPUT} maxworm=240000; {$else INDOUTPUT} maxworm=2000000; {$endif INDOUTPUT} maxwormmargin=maxworm-10000;
      csvheader='  time, av-gen, nworms,  neggs,nlarvae,nadults,ndauers, nbags,      mworms,       meggs,     mlarvae,     madults,     mdauers,       mbags,   Food,  Starve,  LarvFail,  Mature,     Bag,   DieOA,      ...Census of larvae by size... ,,,,,'+
                ' OA@D,    AA@D,    LA@D,    AM@D,    LM@D, TotPred,  Pred-E,  Pred-L,  Pred-D,  Pred-A,  Pred-B,     A@M,    T@A, EggsPerA, RepSpan,' +
                '   FoodL,   FoodA,    XrnL,    XrnA, NewEggs,  Mature,  Awaken, Brt-Dth';
      headerstr='  time|avg gen| nworms|  neggs|nlarvae|nadults|ndauers|  nbags|   Food  Starve:Mature     Bag:DieOA';
{
fi	1	food->in
fp	2	food->predation FP
fx	3	food->expiration FX
el	4	egg->larva
ep	5	egg->predation
fl	6	food->larva
lm	7	larva->maintenance
la	8	larva->adult
ld	9	larva->dauer
lp	10	larva->predation
ls	11	larva->starvation
dt	12	dauer->attrition
dl	13	dauer->larva
dp	14	dauer->predation
fa	15	food->adult
am	16	adult->maintenance
aa	17	adult->aging death
ab	18	adult->bag
ae	19	adult->egg
ap	20	adult->predation
bb	21	bag->burst
bp	22	bag->predation
bd	23	bag->dauer
bw	24	bag->waste
lc	25	larva->CoL
lz	26	larva->efficiency of anabolism
ac	27	adult->CoL
az	28	adult->efficiency of anabolism


Flux Identities

Food:   delta Foodmass = FI – FP – FX – FL – FA
Eggs:   delta Eggmass = AE – EL – EP
        delta Eggcount = AE - EL - EP
Larvae: delta Larvamass = EL4 + FL6 + DL13 – LM7 – LD9 – LA8 – LP10 - LS11
        delta Larvacount = EL – LP – LA – LS
Dauers: delta Dauermass = LD + BD – DL – DP - DT
        delta Dauercount = LD + BD – DL - DP - DT
Adults: delta Adultmass = LA8 + AE15 – AM16 – AA17 – AB18 – AE19 - AP20
        delta Adultcount = LA – AP – AA - AB
Bags:   delta Bagmass =AB – BP – BB – BD
        delta Bagcount = AB – BP – BB
}

      {ratehdr='t,FI,_,FP,_,FX,_,EL,nEL,EP,nEP,EL,nEL,'+
              'FL,_,LM,_,LA,nLA,LD,nLD,LP,nLP,LS,nLS,DT,nDT,DL,nDL,DP,nDP,FA,_,AM,_,AA,nAA,AB,nAB,AE,nAE,AP,nAP,BB,nBB,BP,nBP,BD,nBD,BW';}
      ratehdr='t,FI,FP,FX,EL,nEL,EP,nEP,FL,LM,LA,nLA,LD,nLD,LP,nLP,LS,nLS,DT,nDT,DL,nDL,DP,nDP,FA,AM,AA,nAA,AB,nAB,AE,nAE,AP,nAP,'+
              'BB,nBB,BP,nBP,BD,nBD,BW,LC,LZ,AC,AZ,'+
              'dF,FI-FP-FX-FL-FA,dE,AE-EL-EP,dnE,nAE-nEL-nEP,dL,EL+FL+DL-LD-LM-LP-LA-LS,dnL,EL+DL-LD-LP-LA-LS,dD,LD+BD-DL-DP-DT,dnD,LD+BD-DL-DP-DT,dA,'+
              'LA+FA-AM-AE-AP-AA-AB,dnA,LA-AP-AA-AB,dB,AB-BP-BB-BD,dnB,AB-BP-BB,Food discrep,Egg discrep,,Larva discrep,,Dauer discrep,,Adult discrep,,Bag discrep,,Waste Discrepancies';
      {old ratehdr='t,cr1,cr2,cr3,cr4,cr5,cr6,cr7,cr8,cr9,cr10,cr11,cr12,cr13,cr14,cr15,avgcr1,avgcr2,avgcr3,avgcr4,avgcr5,avgcr6,avgcr7,avgcr8,avgcr9,avgcr10,avgcr11,avgcr12,avgcr13,avgcr14,avgcr15,' +
              'mr1,mr2,mr3,mr4,mr5,mr6,mr7,mr8,mr9,mr10,mr11,mr12,mr13,mr14,mr15,mr16,mr17,mr18,mr19,mr20,mr21,avgmr1,avgmr2,avgmr3,avgmr4,avgmr5,avgmr6,avgmr7,avgmr8,avmgr9,avgmr10,avgmr11,avgmr12,avgmr13,avgmr14,avgmr15,avgmr16,avgmr17avgmr18,avgmr19,avgmr20,avgmr21';}
      rateswithcounts:set of byte=[4,5,8,9,10,11,12,13,14,17..23];
      tempfilename='iamatemporaryfile.csv'; feedfilename:string='FeedSched.inp';   fileincludesculling:boolean=false;
{$ifdef PHEROMONE} woutput:boolean=false; {$else PHEROMONE} woutput:boolean=true; {$endif PHEROMONE}

const hatchage:integer=5;  burstage:integer=10; lifespan:integer=48; {in units of 1/8 day, or 3 hours} fertspan:double=1; {appspan:double=39.2;} {geometric mean betw lifespan and fertspan}  f1=0.6; {const in Gaussian fertility curve}  GompertzN:integer=5;  {controls shape of survival curve, higher number=squarer curve}
      tummult=0.85;
      dietmult:double=0.0;   {adjustment to lifespan based on food consumption}          done:boolean=false;
      {stillbirth:double=0.8;}  minmaturityage:integer=20; maxmaturityage:integer=28;
      foodunit=100; {nanograms - about a portion for an adult}
      CoL:double=0.025; {metabolic cost of living = 2.5% of lean mass per time step}
      dauerfoodthresh:double=250000;  {Dauering}
      maxdaueringrate:double=2;  {number of time periods it takes on average to dauer when there is zero food}   
      dauerawakenprob:double=0.6;
      attrition_time:double=1000; {time scale on which dauers will die}
      larva_app_c0:double=0.196;     larva_app_c1:double=121;           {for calculating appetite function of larvae}
      adult_app_c0:double=0.25;    adult_app_c1:double=0.13508;  adult_food_sat:double=0.460728*5*1E6; {removed c2 and c3 7/15 adult_app_c2:double=0.2;    adult_app_c3:double=15;  {in appetite function for adults}
      minmaturemass:double=1200;  eggmassconst:double=20; newlarvamass:double=65;  dauermassconst:double=400; {will be recalculated as mean proportion betw adult and egg}
      leanmass2eggsproportion:double=0.5; {If adults mature large and are then exposed to low food, they can turn their stored mass into eggs, reducing mass down to a minimum.  This minimum is a weighted average of the mass on attaining maturity and the min lean mass at maturity.  The higher this constant is, the more eggs the stressed adult can lay, turning its lean mass into eggs.}
      fulladultmass:double=3200;
      recipadultmass:double=0.00000000142857;  bagthreshold:double=7000; {Food less than this amount triggers bagging}
      fertconst:double=0.2;       {adults: this much stored energy turned to eggs in each time step}
      anaboliconst:double=0.63212;   {larvae: this much stored energy turned to lean mass in each time step}
      efficiency:double=0.75;  {efficiency of converting food to leanmass and to eggs, used to be separate, now just one}
      pumplimit:double=0.001;           {If food is scarce, what portion of total food in the test tube can a worm eat in a period?}
      initeggs:integer=2000;
      initdauers:integer=500;
      initlarvae:integer=0;
      {initlarvae:integer=2000;
      initwormbags:integer=200;
      initadults:integer=500;}
      foodadded:double=100000;
      foodexpire:integer=10000;
      foodattrition:double=0.99;
      feedinginterval:integer=16;
      cullinginterval:integer=16;
      cullportion:double=0.2;          {this proportion of worms is "predated" with every feeding}
      foodcullportion:double=0.0;         {2016Oct at Kerry's request, predation proportion is individualized by 5 stages + food}
      eggcullportion:double=0.0;
      larvacullportion:double=0.0;
      dauercullportion:double=0.0;
      adultcullportion:double=0.0;
      bagcullportion:double=0.0;
      stoptime:integer=3000;
      feedingdisparity:double=0;
      escapestarvationprob:double=0;
      headerepeat:integer=50;
      arraysiz=30000;
      lookback=0.94;  oneminuslookback=1-lookback;  {constants for averaging rates}
      randomfeeding:boolean=false; maxterval:word=16;  minterval:word=16;
      scheduledfeeding:boolean=false;
      randomculling:boolean=false; maxcerval:word=16;  mincerval:word=16;
      feedandculltogether:boolean=true;
      stopper:boolean=false;
      synchronized:boolean=false;
	  bm=0.000129574732510846;  bn=0.000248858348410121*5E6;  {Kr=1.83115642756057; Ks=-0.244862869993007*5E6;  {4 parameters apperaing in the logistic eqn for worm growth, derived in Modeling_Growth.xlsx.  Derived from these are b and K, which also depend on food.
	  The 5E6 is because food in the spreadsheet is in mg/ml, while food in this program is in nanograms per 5ml.  a depends on b and new-hatched larva mass.}
	  {K=Kr+Ks/food; has been replaced by K=Kr*tanh(Ks*food)} Kr=1.78027908103543/8; {1/8 converts from days to pds}  Ks=2.13217438144031/5E6;

type  sixdoubles=array[1..6] of double;
const
      fd:sixdoubles=(5E6*0.07,5E6*0.13,5E6*0.25,5E6*0.50,5E6*1.0,5E6*4.0);
      xp:sixdoubles=(1.319189633,2.201468416,3.484973388,3.640140648,4.725008562,4.156785526);
      agescale:sixdoubles=(3.481473101,2.024212254,1.111015285,0.897214743,0.739552921,0.749568986);
      fertmult:sixdoubles=(3.272545277/8,4.817111781/8,7.326346126/8,15.73194284/8,8.617245216/8,16.06551326/8);

      {fdata:array[0..3] of sixdoubles=((5E6*0.07,5E6*0.13,5E6*0.5,5E6*1.0),(4.156785526,4.725008562,3.484973388,1.319189633),(0.749568986*8,0.739552921*8,1.111015285*8,3.481473101*8),(16.06551326,8.617245216,7.326346126,3.272545277));}

type doppelganger=^double;
     intpointer=^integer;
     intarray=array[1..arraysiz] of integer;

{$ifdef INDOUTPUT}
     onewormonetimestep=record
       foodavail                 :single;
       leanm,{fatm,}
       foodeaten                 :word;
       eggslaid                  :byte;
       state                     :char; {L=larva; D=dauer; W=wakeup; M=mature; B=bag; R=burst; P=predated; +=death}
       end;
     timesteps=array[0..maxlife] of onewormonetimestep;
     tsptr=^timesteps;
{$endif INDOUTPUT}

     worm=object
       s         :(egg, larva, dauer, adult, wormbag, dead);
       birthday  :single;  {when egg laid} {also useful as pointer to the lifehistory file}
       tatlastegg,tatfirstegg :single; {when did I lay my first egg? my last?}
       ttime     :single;  {“transition time” when did the worm enter its current state}
       minreducedmass :single;  {lean mass on attaining maturity}
       leanmass, {fatmass,} eggmass  :single; {total body energy is the sum of lean mass + reproductive tissue + stored energy}
       lastfood,           {memory of food available in last time step}
       avgadultfood,  {memory of food history during adult stage only}
       totfood     :single;  {total food consumed over a lifetime, as a proportion of the max that the worm could have eaten if limitless food were available}
       enterdauer_t  :single; {time at which larva entered dauer state; 0 means it was never a dauer}
       cum_eggs  :word;
       gen       :word;
       age       :smallint;
       frombag   :boolean;
{$ifdef VARIABILITY}
       ipumplimit,iminmaturemass,imaturityage,ilifespan   :single;
{$endif VARIABILITY}
{$ifdef INDOUTPUT}
        modeofdeath :char; {P=predation, B=bag,  A=oldage,  S=starvation,  T=attrition}
        matur       :boolean;
        lh          :tsptr;
       procedure RecordLHState;
       function stageLetter:char;
       function AdultAge:smallint;
{$endif INDOUTPUT}
       function Appetite(var fertrq:double):double;
       function Stage:string;
       {function age :smallint; replaced with variable}
       procedure Eat;
       procedure Hatch;
       function rage:double;
       function stagespecificcull:double;
       function Gompertz:double; {age-adjusted mortality}
       function AgeAdjustedFertility(setback :double):double;
       {function AgeAdjustedAppetite:double;}
       procedure Mature;
       function Siz:double; {length as a function of lean mass}
       function Attrition:boolean; {dauer}
       procedure Awaken; {dauer}
       procedure Reproduce;
       procedure NewWorm(gg :integer);
       procedure Alchemy; {energy stored as fat is turned to lean mass in a larva, into egg mass in an adult}
       function Bag:boolean;
       procedure Burst;
       procedure Die;
       function Punt:boolean; {larva -> dauer}
       procedure Starve; {larva death}
       procedure Predate;
       procedure Eliminate;
       procedure SenescentMortality;
{$ifdef INDOUTPUT}
       procedure WriteOneLH;
       procedure OneSampleFileLine(adultspan :smallint);
{$endif INDOUTPUT}
       end; {Worm object}

     loop=record
       case isdouble:boolean of
         true :
                (max  :byte;
                 linenum :byte;
                 valu :array[1..16] of double;
                 dbladr  :doppelganger);
        false :
                (dummy :byte;
                 dummy2 :byte;
                 dummy3 :array[1..16] of double;
                 intadr  :intpointer);
       end;

var fin,fout,ratesfile,samplefile,feedsched     :text;
    samplen,samplecount,nextfeeding,
    PRcount,peakPRcount,refPRcount,
    allworms,matureworms,oldageworms            :integer;
    sampletype                                  :char; {A=all, M=reached maturity, O=died of old age, E=Early, F=Fertile}
    pop                                         :intarray;
    lastpd                                      :word;
    t,lastfeeding,lastculling                   :double; {in units of 3 hrs, but incremented continuously}
    foutname,barefoutname,finname,ratefname     :string;
    w                                           :array[1..maxworm] of worm;
    census                                      :array[1..5] of integer;
{$ifdef VARIABILITY}
    indvariation                                :double;
{$endif VARIABILITY}
    food,oldfood,rememberedfood                 :double; {extensive for the whole tube, measured in nanograms}
    appsum                                      :double; {total pumping capacities of all worms and larvae, used for pro rata distribution of food}
    nworms,neggs,nlarvae,nadults,ndauers,nbags,
    oldnworms,oldneggs,oldnlarvae,oldnadults,
    oldndauers,oldnbags                         :integer;
    mworms,{meggs,}mlarvae,madults,mdauers,mbags,
    oldmworms,{oldmeggs,}oldmlarvae,oldmadults,
    oldmdauers,oldmbags                         :double;
    starvect,fadeawayct,maturect,bagct          :integer;
    overflow,sampleon                           :boolean;
    LL                                          :array[1..4] of loop;
    olda,adulta,larvaa,adultm,larvam            :double; {age and mass at death}
    oldct,adultct,larvact                       :integer; {death counts}
    predct,pred_e,pred_l,pred_d,pred_a,pred_b   :integer; {number predated: egg,larva,dauer,adult,bag}
    AatM,TasA,EggsPerA,RepSpan                  :double;  {average age at maturity; average time spent as adult; average eggs laid among adults that died this period; fert span = average time between first and last egg laid by adults that died this period}
    NewEggs,Awakenct,BirthsMinusDeaths          :integer;
    countr                                      :array[0..28] of integer;{rate output, # of worms} {mass and count rates are now indexed to same number, 2016Aug}
    {avcountr                                    :array[0..24] of double; {rate output, # of worms}
    massr{,avmassr}                             :array[0..28] of double; {rate output, biomass}
    mean,stdev                                  :double;
    tau, A                                      :double; {constants in Gompertz curve}
  {$ifdef INDOUTPUT}
    invdirectory                                :string;
{$endif INDOUTPUT}
{$ifdef PHEROMONE}
    PheroFile                                   :text;
    dauercount                                  :integer;
    adjustedfoodthreshold,
    lambda,mu,persistencetime,
    cumagesum,avgcumagesum                      :double;
{$endif PHEROMONE}
{$ifdef NDEBUG}
function randomn(r :integer):integer;
         const last :double=0.5;
         begin
         result:=round(r*last);
         last:=1.3*result/r;
         if (last>1) then last:=last-1;
         end;

function randomr:double;
         const last :double=0.5;
         begin
         result:=last*1.3;
         if (result>1) then result:=result-1;
         end;

function FastPoissonDist(r :double):integer;
         begin
         result:=random(round(int(r))); if (random<r-int(r)) then inc(result);
         end;
{$else}
function randomn(r :integer):integer;
         begin
         result:=random(r);
         end;

function randomr:double;
         begin
         result:=random;
         end;
{$endif}

procedure LarvaCensus;
          const crit:array[1..4] of double=(6.9314484316, 24.0224886796, 83.2553207402, 288.5399811814);
          var i,j: integer;
          begin
          for j:=1 to 5 do census[j]:=0;
          for i:=1 to nworms do with w[i] do if (s=larva) then begin
             j:=0;
             repeat inc(j) until (j=5) or (crit[j]>leanmass);
             inc(census[j]);
          end;end;

procedure ReadFromSchedFile;
          begin
          if (eof(feedsched)) then begin nextfeeding:=maxint; exit; end;
          if (fileincludesculling) then readln(feedsched,nextfeeding,foodadded,cullportion)
          else readln(feedsched,nextfeeding,foodadded);
          end;

procedure ReadInputFile;
          var loopi,line :byte;
              cullx      :double;

  procedure ReadBoolean(var b :boolean);
            var s  :string;
                ch :char;
            begin
            s:='';
            repeat read(fin,ch); s:=s+upcase(ch); until eof(fin) or (ch in ['A'..'z']);
            repeat read(fin,ch); s:=s+upcase(ch); until eof(fin) or (pos(s,'TRUE')>0) or (pos(s,'FALSE')>0);
            if (pos(s,'FALSE')>0) then b:=false
            else if (pos(s,'TRUE')>0) then b:=true;
            end;

  procedure ReadFileName;
            var ch :char;
            begin
            foutname:='';
            repeat read(fin,ch); until eof(fin) or (ch='''');
            repeat read(fin,ch); foutname:=foutname+upcase(ch); until eof(fin) or (ch='''');
            SetLength(foutname,pred(length(foutname)));
            barefoutname:=foutname;
            end;

  procedure ReadMultipleDoubles(v :doppelganger);
            var ch :char;
            begin
            read(fin,v^);
            repeat read(fin,ch) until (not (ch in whitespace));
            if (ch<>'&') then exit;
            if (loopi=4) then begin writeln('Number of loop variables is limited to 4.'); write('Press <ENTER>, then change input file.'); readln; halt; end;
            inc(loopi);
            with LL[loopi] do begin
               isdouble:=true;
               linenum:=line;
               dbladr:=v; valu[1]:=v^;
               repeat
                 inc(max); read(fin,valu[max]); repeat read(fin,ch) until (not (ch in whitespace));
               until (ch<>'&') or (max=16);
               end;
            end;

  procedure ReadMultipleInts(v :intpointer);
            var ch :char;
            begin
            read(fin,v^);
            repeat read(fin,ch) until (not (ch in whitespace));
            if (ch<>'&') then exit;
            if (loopi=4) then begin writeln('Number of loop variables is limited to 4.'); write('Press <ENTER>, then change input file.'); readln; halt; end;
            inc(loopi);
            with LL[loopi] do begin
               isdouble:=false;
               linenum:=line;
               intadr:=v; valu[1]:=v^;
               repeat
                 inc(max); read(fin,valu[max]); repeat read(fin,ch) until (not (ch in whitespace));
               until (ch<>'&') or (max=16);
               end;
            end;

  procedure ReadSampleInstructions;
            var ch :char;
            begin
            read(fin,samplen);
            repeat read(fin,ch) until (ch>' ');
            if (ch in ['E','A','O','M','F']) then sampletype:=ch else sampletype:='A';
            end;

  procedure ReadParameter;
            var s  :string;
                ta :smallint;
                ch :char;
            begin
            inc(line);
            s:='';
            repeat read(fin,ch); s:=s+ch; until eof(fin) or (ch='=') or (ch='{');
            if (ch='{') then begin readln(fin,s); exit; end; {Ignore lines that begin with bracket}
            TrimToUpper(s);
            if (pos('OUTPUT FILE NAME',s)>0) then ReadFileName
            else if (pos('INITIAL NUMBER OF EGGS',s)>0) then ReadMultipleInts(@initeggs)
            else if (pos('INITIAL NUMBER OF DAUERS',s)>0) then ReadMultipleInts(@initdauers)
            else if (pos('INITIAL NUMBER OF LARVAE',s)>0) then ReadMultipleInts(@initlarvae)
            else if (pos('AMOUNT OF FOOD ADDED',s)>0) then ReadMultipleDoubles(@foodadded)
            else if (pos('TIME BETWEEN FEEDINGS',s)>0) then ReadMultipleInts(@feedinginterval)  {7223 is a special signal, meaning Mon Wed Fri feeding}
            else if (pos('TIME BETWEEN CULLINGS',s)>0) then ReadMultipleInts(@cullinginterval)  {7223 is a special signal, meaning Mon Wed Fri feeding}
            else if (pos('DAUER FOOD THRESHOLD',s)>0) then ReadMultipleDoubles(@dauerfoodthresh)
            else if (pos('INVERSE DAUERING RATE',s)>0) then ReadMultipleDoubles(@maxdaueringrate)
            else if (pos('DAUER ATTRITION',s)>0) then ReadMultipleDoubles(@attrition_time)
            else if (pos('CONST FOR TOTAL FOOD IN DEATH PROB',s)>0) then ReadMultipleDoubles(@dietmult)
            else if (pos('TIME FOR EGG TO HATCH',s)>0) then ReadMultipleInts(@hatchage)
            else if (pos('COST OF LIVING',s)>0) then ReadMultipleDoubles(@CoL)
            else if (pos('MIN LEAN MASS AT MATURITY',s)>0) then ReadMultipleDoubles(@minmaturemass)
            else if (pos('MAX LEAN MASS AT MATURITY',s)>0) then ReadMultipleDoubles(@fulladultmass)
            else if (pos('PROPORTION CONST FOR MIN MASS OF STRESSED ADULTS',s)>0) then ReadMultipleDoubles(@leanmass2eggsproportion)
            {else if (pos('PROPORTION LEAN MASS GROWTH DURING ADULT PHASE',s)>0) then ReadMultipleDoubles(@adultmassgrowthpct)}
            else if (pos('MASS OF ONE EGG',s)>0) then ReadMultipleDoubles(@eggmassconst)
            else if (pos('MASS OF NEWLY HATCHED LARVA',s)>0) then ReadMultipleDoubles(@newlarvamass)
            {else if (pos('LEAN MASS OF DAUER',s)>0) then ReadMultipleDoubles(@dauermassconst) fixed at mean proportion between egg and adult}
            {else if (pos('STARVATION THRESHOLD',s)>0) then ReadMultipleDoubles(@larvastarvethresh)  fixed to zero}
            else if (pos('PROBABILITY OF DAUER AWAKENING',s)>0) then ReadMultipleDoubles(@dauerawakenprob)
            else if (pos('MIN FOOD REQUIRED TO AVOID BAGGING',s)>0)  then ReadMultipleDoubles(@bagthreshold)
            else if (pos('LARVA ANABOLIC CONST',s)>0) then ReadMultipleDoubles(@anaboliconst)
            else if (pos('METABOLIC EFFICIENCY',s)>0) then ReadMultipleDoubles(@efficiency)
            else if (pos('DEGRADATION TIME FOR FOOD',s)>0) then ReadMultipleInts(@foodexpire)
            else if (pos('MIN AGE AT MATURITY',s)>0) then ReadMultipleInts(@minmaturityage)
            else if (pos('EXPIRY DATE FOR LARVAE',s)>0) or (pos('MAX AGE AT MATURITY',s)>0) or (pos('TIME TARGET FOR LARVA MATURITY',s)>0) then ReadMultipleInts(@maxmaturityage)
            else if (pos('ADULT REPRODUCTIVE CONST',s)>0) then ReadMultipleDoubles(@fertconst)
            else if (pos('LARVA APPETITE CONST #0',s)>0) then ReadMultipleDoubles(@larva_app_c0)
            else if (pos('LARVA APPETITE CONST #1',s)>0) then ReadMultipleDoubles(@larva_app_c1)
            else if (pos('ADULT APPETITE CONST #0',s)>0) then ReadMultipleDoubles(@adult_app_c0)
            else if (pos('ADULT APPETITE CONST #1',s)>0) then ReadMultipleDoubles(@adult_app_c1)
            else if (pos('ADULT FOOD SATURATION',s)>0) then ReadMultipleDoubles(@adult_food_sat)
            {else if (pos('ADULT FOOD FILTRATION CONST',s)>0) then ReadMultipleDoubles(@pumplimit) {same for adults and larvae}
            else if (pos('LARVA FOOD FILTRATION CONST',s)>0) then ReadMultipleDoubles(@pumplimit) {applies only to larvae; adults are governed by "Adult food saturation"}
            else if (pos('LIFE SPAN',s)>0) then ReadMultipleInts(@lifespan)
            else if (pos('SHAPE OF SURVIVAL CURVE',s)>0) then ReadMultipleInts(@GompertzN)
            else if (pos('FERTILITY SPAN',s)>0) then ReadMultipleDoubles(@fertspan)
            else if (pos('TIME FROM BAGGING TO BURSTING',s)>0) then ReadMultipleInts(@burstage)
            else if (pos('OVERALL CULLING PROPORTION (PREDATION)',s)>0) then ReadMultipleDoubles(@cullportion)
            else if (pos('FOOD CULLING PROPORTION (PREDATION)',s)>0) then ReadMultipleDoubles(@foodcullportion)
            else if (pos('EGG CULLING PROPORTION (PREDATION)',s)>0) then ReadMultipleDoubles(@eggcullportion)
            else if (pos('LARVA CULLING PROPORTION (PREDATION)',s)>0) then ReadMultipleDoubles(@larvacullportion)
            else if (pos('DAUER CULLING PROPORTION (PREDATION)',s)>0) then ReadMultipleDoubles(@dauercullportion)
            else if (pos('ADULT CULLING PROPORTION (PREDATION)',s)>0) then ReadMultipleDoubles(@adultcullportion)
            else if (pos('BAG CULLING PROPORTION (PREDATION)',s)>0) then ReadMultipleDoubles(@bagcullportion)
            else if (pos('STOP AFTER HOW MANY TIME STEPS',s)>0) then read(fin,stoptime)
            else if (pos('REPEAT HEADER EVERY N TIME STEPS',s)>0) then read(fin,headerepeat)
            else if (pos('INDIVIDUAL WORM LIFE HISTORY SAMPLE',s)>0) then ReadSampleInstructions
            else if (pos('FOOD TRIAGE',s)>0) then ReadMultipleDoubles(@feedingdisparity)
            else if (pos('LARVA PROBABILITY TO CHEAT STARVATION',s)>0) then ReadMultipleDoubles(@escapestarvationprob)
{$ifdef VARIABILITY}
            else if (pos('INDIVIDUAL VARIATION IN PUMPING AND ADULT SIZE',s)>0) then ReadMultipleDoubles(@indvariation)
{$else VARIABILITY}
            else if (pos('INDIVIDUAL VARIATION IN PUMPING AND ADULT SIZE',s)>0) then 
{$endif VARIABILITY}

            else if (pos('SCREEN OUTPUT',s)>0) then ReadBoolean(woutput)
            else if (pos('SYNCHRONIZED',s)>0) then ReadBoolean(synchronized)
            else if (s>'') then begin writeln('Unrecognized parameter: ',s); readln; end;
           {$I-} readln(fin,s); {$I+} ta:=ioresult; {throw away}
            end;

          begin {ReadInputFile}
{$ifdef PHEROMONE} finname:='Pheromone.inp'; {$else} finname:='WormNRG.inp'; {$endif PHEROMONE} {finname:='merge.inp';}
          foutname:='WormNRG.csv'; ratefname:='WRates.csv';
          assign(fin,finname); {$I-} reset(fin); {$I+}
          if (ioresult>0) then begin
            write('File not found: "',finname);
            writeln('".  Using default parameters.');
            exit;
            end;
          for loopi:=1 to 4 do with LL[loopi] do begin max:=1; linenum:=0; end;
          loopi:=0; line:=0;
          repeat ReadParameter until eof(fin);
          close(fin);
          maxdaueringrate:=one/maxdaueringrate; {number specified in inp file is # of periods} 
          dauermassconst:=sqrt(newlarvamass*minmaturemass);
          bagthreshold:=2.0*bagthreshold; {because it will be compared to the sum of two time steps}
          tau:=tummult*lifespan/GompertzN;   A:=lifespan*(exp(GompertzN)-one);    {constants in Gompertz curve}
          if (foodexpire=0) then foodattrition:=1 else foodattrition:=exp(-1/foodexpire);
          {no more appspan, 9/17 appspan:=sqrt(fertspan*lifespan); {appetite falls at an age intermediate between fertility and mortality}
          if (headerepeat=0) then headerepeat:=maxint;
          feedandculltogether:= (cullinginterval=feedinginterval);
          if (feedinginterval=-1) then begin
            scheduledfeeding:=true;
            {$I-}assign(feedsched,feedfilename); reset(feedsched);{$I+}
            if (ioresult>0) then begin write('"FeedSched.inp" file not found. Press <Enter> to quit.'); readln; halt; end;
            read(feedsched,nextfeeding,foodadded);
            readln(feedsched,cullx); if (cullx<1) then begin fileincludesculling:=true; feedandculltogether:=true; end;
            reset(feedsched);
            ReadFromSchedFile;
            end;
          if (feedinginterval>100) and (feedinginterval<>7322) and (feedinginterval<>7223) then begin
            randomfeeding:=true; minterval:=feedinginterval div 100;  maxterval:=feedinginterval mod 100;
            if (minterval>maxterval) then {swap} begin feedinginterval:=maxterval; maxterval:=minterval; minterval:=feedinginterval; end;
            feedinginterval:=maxterval+random(succ(maxterval-minterval));
            end;
          if (not feedandculltogether) then begin
             if (cullinginterval=0) then cullinginterval:=maxint
             else if (cullinginterval>100) and (cullinginterval<>7322) and (cullinginterval<>7223) then begin
               randomculling:=true; mincerval:=cullinginterval div 100;  maxcerval:=cullinginterval mod 100;
               if (mincerval>maxcerval) then {swap} begin cullinginterval:=maxcerval; maxcerval:=mincerval; mincerval:=cullinginterval; end;
               cullinginterval:=maxcerval+random(succ(maxcerval-mincerval));
               end;
             end;
          end; {ReadInputFile}

{total food=200000*4000/16= 5E7; total worms=136528;  mature worms=9466; died of old age=3522}

procedure OpenOutputFile(i1,i2,i3,i4 :byte);
          var which :array[1..4] of byte;

   procedure OutputEditedLine(i :byte; var copyto :text);
             var j :byte;
                ch :char;
                ws :string;
             begin
             repeat read(fin,ch); write(copyto,ch); until (ch='=');
             if (which[i]>1) then for j:=1 to which[i]-1 do repeat read(fin,ch) until (ch='&');
             repeat read(fin,ch); write(copyto,ch); until (ch in numchars);
             repeat read(fin,ch); write(copyto,ch); until (not (ch in numchars));
             readln(fin,ws);   if (pos('{',ws)>0) then writeln(copyto,copy(ws,pos('{',ws),length(ws))) else writeln(copyto);
             end;

   procedure CopyInputToOutputFile(var copyto :text);
             var ws :string;
                 line,i :byte;
             begin
             reset(fin); line:=0;
             writeln(copyto,crlf,vxstr + DateTimeToStr(date+time){$ifdef VARIABILITY}+' VARIABILITY' {$endif});
             if (synchronized) then writeln(copyto,'SYNCHRONIZED');
             while not eof(fin) do begin
               inc(line);
               i:=0; repeat inc(i) until (i>4) or (LL[i].linenum=line);
               if (i>4) then begin readln(fin,ws); writeln(copyto,ws); end
               else OutputEditedLine(i,copyto);
               end;
             close(fin);
             end;

          begin
          which[1]:=i1; which[2]:=i2; which[3]:=i3; which[4]:=i4;
          {assign(tempfile,'WormNRG.tmp');  rewrite(tempfile); close(tempfile);}
          {foutname:='WormNRG.csv';}
          assign(fout,foutname); if (fileexists(foutname)) then append(fout) else rewrite(fout);
          if true then CopyInputToOutputFile(fout)
          else begin
            writeln(fout,'Init worm ct=',initeggs+initdauers+initlarvae,'; Time betw feedings=',feedinginterval,'; Time betw cullings=',cullinginterval,'; Food added=',foodadded,
                    '; lifespan=',lifespan);
            writeln(fout,'Appetite constants for larva=',larva_app_c0:5:2,larva_app_c1:5:2,'  for adults=',adult_app_c0:5:2,adult_app_c1:5:2);
            writeln(fout,'Time for bag to burst=',burstage);
            writeln(fout,' Max proportion of food a worm can gather=',FtoaNoE(pumplimit,5,3),'; Culling pct=',round(100*cullportion),'%');
            end;
          close(fout);
          assign(fout,tempfilename); rewrite(fout);
          if (upcase(foutname[length(foutname)])='V') then writeln(fout,CSVheader) else writeln(fout,headerstr);
          close(fout);
          assign(ratesfile,ratefname); {$I-} append(ratesfile); {$I+} if (ioresult>0) then rewrite(ratesfile);
          CopyInputToOutputFile(ratesfile);
          writeln(ratesfile,ratehdr);
          close(ratesfile);
          end;

procedure InitializeCountsToZero;
          var i :byte;
          begin
          predct:=0; pred_e:=0; pred_l:=0; pred_d:=0; pred_a:=0; pred_b:=0;
          AatM:=0; TasA:=0; EggsPerA:=0; RepSpan:=0;
          starvect:=0; fadeawayct:=0; maturect:=0; bagct:=0; oldct:=0;
          olda:=0; adulta:=0; larvaa:=0; adultm:=0; larvam:=0; oldct:=0; adultct:=0; larvact:=0;
          NewEggs:=0; Awakenct:=0; BirthsMinusDeaths:=0;
          for i:=1 to 28 do begin countr[i]:=0; massr[i]:=0; end;
          end;

procedure ComposeSampleFileHeader;
          var i :word;
              ws :string;
          begin
          reset(fin);
          writeln(samplefile,vxstr + DateTimeToStr(date+time));
          while not eof(fin) do begin
            readln(fin,ws); writeln(samplefile,ws);
            end;
          close(fin);
          write(samplefile,'BDay, Age@Death, Mode_death');
          for i:=0 to maxlife do write(samplefile,', _',i,'_, Lean_mass,Fat_mass,Food_avail,Food_eaten,Eggs_laid');
          writeln(samplefile);
          end;

procedure Init;
          var i         :integer;
          begin
          {$ifndef NDEBUG} randomize; {$endif NDEBUG}
          allworms:=0; matureworms:=0; oldageworms:=0; peakPRcount:=0; refPRcount:=maxint;
          dauermassconst:=sqrt(minmaturemass*eggmassconst);  recipadultmass:=1/minmaturemass;
          sampleon:=(sampletype='E');
          InitializeCountsToZero;
          {for i:=1 to 24 do avcountr[i]:=0;
          for i:=1 to 24 do avmassr[i]:=0;  {averages are initialized only once}
          t:=0; lastpd:=0;
          nworms:=initeggs+initdauers+initlarvae; food:=foodadded; oldfood:=food; lastfeeding:=0; lastculling:=0;
          oldnworms:=0; oldneggs:=initeggs; oldnlarvae:=0; oldndauers:=initdauers; oldnadults:=0; oldnbags:=0;
          oldmworms:=0; oldmlarvae:=0; oldmdauers:=0; oldmadults:=0; oldmbags:=0;
          {for i:=1 to maxworm do w[i].lh:=new(tsptr);}
{$ifdef INDOUTPUT}
          for i:=1 to 30000 do w[i].lh:=new(tsptr);
          for i:=30001 to maxworm do w[i].lh:=nil;  {shouldn't be necessary}
{$endif INDOUTPUT}
          for i:=1 to initeggs do with w[i] do begin
            gen:=0;
            s:=egg;
            cum_eggs:=0;
            tatfirstegg:=0; tatlastegg:=0;
            birthday:=-random*hatchage;
            if (synchronized) then birthday:=-hatchage-1+tenth*(random-random); {so they'll all hatch immediately}
            ttime:=birthday;
            enterdauer_t:=0;
            totfood:=0;
            avgadultfood:=0;
            leanmass:=eggmassconst;
            {fatmass:=0;}
{$ifdef VARIABILITY}
            ipumplimit:=exp(GaussRand*indvariation)*pumplimit;
            iminmaturemass:=exp(GaussRand*indvariation)*maturemassmin;
            imaturityage:=exp(GaussRand*indvariation)*minmaturityage;
            ilifespan:=exp(GaussRand*indvariation)*lifespan;
{$endif VARIABILITY}
            end;
         for i:=succ(initeggs) to (initeggs+initdauers) do with w[i] do begin
            gen:=0;
            s:=dauer;
            cum_eggs:=0;
            tatfirstegg:=0; tatlastegg:=0;
            birthday:=-(40*random)*hatchage;
            ttime:=birthday+hatchage+random(10);
            enterdauer_t:=t;
            totfood:=-ttime*random*foodunit;
            avgadultfood:=0;
            leanmass:=dauermassconst;
            {fatmass:=half*dauermassconst;}
            eggmass:=0;
            oldmdauers:=oldmdauers+leanmass {+fatmass};
{$ifdef VARIABILITY}
            ipumplimit:=exp(GaussRand*indvariation)*pumplimit;
            iminmaturemass:=exp(GaussRand*indvariation)*maturemassmin;
            imaturityage:=exp(GaussRand*indvariation)*minmaturityage;
{$endif VARIABILITY}
            end;
          for i:=succ(initeggs+initdauers) to nworms do with w[i] do begin
            gen:=0;
            s:=larva;
            cum_eggs:=0;
            tatfirstegg:=0; tatlastegg:=0;
            birthday:=-(3*random)*hatchage;
            ttime:=birthday+hatchage+random(3);
            enterdauer_t:=0;
            totfood:=-ttime*random*foodunit;
            avgadultfood:=0;
            leanmass:=dauermassconst;
            {fatmass:=half*dauermassconst;}
            eggmass:=0;
            oldmlarvae:=oldmlarvae+leanmass{+fatmass};
{$ifdef VARIABILITY}
            ipumplimit:=exp(GaussRand*indvariation)*pumplimit;
            iminmaturemass:=exp(GaussRand*indvariation)*maturemassmin;
            imaturityage:=exp(GaussRand*indvariation)*minmaturityage;
{$endif VARIABILITY}
            end;
(*
          for i:=succ(initeggs+initlarvae) to initeggs+initlarvae+initadults do with w[i] do begin
            gen:=0;
            s:=adult;
            cum_eggs:=10;
            tatfirstegg:=0; tatlastegg:=0;
            birthday:=-(one+random)*hatchage-8;
            ttime:=birthday+hatchage+8;
              if (ttime>0) then begin
                 writeln(birthday, hatchage:8);
                 writeln(ttime);
                 readln;
                 end;
            dauermem:=false;
            totfood:=-ttime*random*foodunit;
            leanmass:=adultmassconst;
            fatmass:=adultmassconst*0.3;
            oldmadults:=oldmadults+leanmass{+fatmass};
            end;
          for i:=succ(initeggs+initlarvae+initadults) to nworms do with w[i] do begin
            gen:=1;
            s:=wormbag;
            birthday:=-(one+random)*hatchage-hatchage;
            tatfirstegg:=0; tatlastegg:=0;
            dauermem:=false;
            ttime:=-random*burstage;
            totfood:=random*foodunit;
            leanmass:=adultmassconst;
            eggmass:=eggmassconst*efficiency*leanmass/dauermassconst; {remember for when it bursts}
            neggs:=round(eggmass/eggmassconst);
            oldmbags:=oldmbags+leanmass{+fatmass}+eggmass;
            end;
*)
          if (samplen>0) then begin
            if (sampletype='E') then samplecount:=0 else samplecount:=maxint; {so we don't start collecting data yet}
            assign(samplefile,'IndividualWorms.csv'); rewrite(samplefile);
            ComposeSampleFileHeader;
            close(samplefile);
            end;
{$ifdef PHEROMONE}
            persistencetime:=5; {1/e-life of hypothetical PRLS pheromone}
            cumagesum:=0; for i:=1 to nworms do with w[i] do
              if (s=adult) then cumagesum:=cumagesum+Rage;
            if (cumagesum=0) then cumagesum:=1;
            avgcumagesum:=cumagesum;
            dauercount:=0;
{$endif PHEROMONE}
          if (nextfeeding=0) and (feedinginterval=-1) then ReadFromSchedFile;
{$ifdef INDOUTPUT}
          DateTimetoString(invdirectory,'yymmdd@hhnn',date+time);
          {$I-} MkDir(invdirectory); {$I+} ioresult;
          invdirectory:=invdirectory+'\W';
{$endif INDOUTPUT}
          end; {Init}

{$ifdef PHEROMONE}
procedure UpdateCumagesum;
          var i      :integer;
              agesum :double;
          begin
          agesum:=0;
          for i:=1 to nworms do with w[i] do
            if (s=adult) then agesum:=agesum+Rage;
          cumagesum:=(agesum+cumagesum*(persistencetime-1)) / persistencetime;
          avgcumagesum:=(cumagesum + avgcumagesum*(t-1)) / t;
          end;
{$endif PHEROMONE}

(*
function worm.age :smallint;
         var r :smallint;
         begin
         r:=floor(t)-floor(birthday);
         if (not frombag) then dec(r,hatchage); {egg hatching time}
         if (r<1) then result:=1
         else result:=r;
         end;
*)
{$ifdef INDOUTPUT}

procedure worm.OneSampleFileLine(adultspan :smallint);
          var i :smallint;
          begin
          write(samplefile,ftoanoE(birthday,0,4),',',adultspan,',',modeofdeath);
          i:=-1;
          repeat
            inc(i);
            with lh^[i] do write(samplefile,',',lh^[i].state,',',ftoanoE(lh^[i].leanm,0,4),',',{ftoanoE(lh^[i].fatm,0,4),',',}round(lh^[i].foodavail),',',ftoanoE(0.1*lh^[i].foodeaten,0,4),',',lh^[i].eggslaid);
          until (lh^[i].state='+') or (lh^[i].state='P') or(i=maxlife);
          writeln(samplefile);
          end;

procedure worm.WriteOneLH;
          var i :smallint;
              invfile :text;
              birthmodeanddauer :string[12];
              costofliving,
              larvalspan,dauerspan,adultspan,fertilespan :byte;
          begin
          if (birthday<100) and (sampletype<>'E') then exit;

          if (frombag) then birthmodeanddauer:='BAG'
          else if (enterdauer_t<>0) then birthmodeanddauer:='EGG to DAUER' else birthmodeanddauer:='EGG';
          if (frombag) then birthmodeanddauer:='BAG';
          larvalspan:=0; dauerspan:=0; adultspan:=0; fertilespan:=0; i:=0;
          repeat
            inc(i);
            case (lh^[i].state) of
               'L' : inc(larvalspan);
               'A' : inc(adultspan);
               'D' : inc(dauerspan);
               'X' : begin write('Missing LH record ',i,' worm #',birthday:7:3);
                     writeln('Birthday=',ftoanoE(birthday,0,4),'; from ',birthmodeanddauer,'; Age at death=',adultspan,'; Mode of death="',modeofdeath,'"');
                     writeln('Larval span=',larvalspan,'; Dauer span=',dauerspan,'; Adultspan=',adultspan,'; Fertile span=',fertilespan,'; Lifetime eggs=',cum_eggs,crlf);
                     OneSampleFileLine(-1);
                     exit;
                     end;
               end;
            if (lh^[i].eggslaid>0) then inc(fertilespan);
          until (lh^[i].state='+') or (lh^[i].state='P') or (i=maxlife);
          if (adultspan>0) then inc(adultspan);
          OneSampleFileLine(adultspan);
          assign(invfile,invdirectory+StringReplace(ftoanoE(birthday,0,3),'.','-',[rfReplaceAll])+'.csv'); rewrite(invfile);
          writeln(invfile,'Birthday=',ftoanoE(birthday,0,4),'; from ',birthmodeanddauer,'; Age at death=',adultspan,'; Mode of death="',modeofdeath,'"');
          writeln(invfile,'Larval span=',larvalspan,'; Dauer span=',dauerspan,'; Adultspan=',adultspan,'; Fertile span=',fertilespan,'; Lifetime eggs=',cum_eggs,crlf);
          writeln(invfile,'age,state,lean-m,delta-lean,delta-fat,food-avail,food-eaten,CoL,#eggs-laid,m-eggs-laid');
          i:=-1; i:=0;
          repeat
            inc(i);
            with lh^[i] do begin
            if (state='L') or (state='A') then costofliving:=round(CoL*leanm) else costofliving:=0;
            writeln(invfile,i,',',state,',',leanm,',',leanm-lh^[pred(i)].leanm,',',{round(fatm),',',round((fatm-lh^[pred(i)].fatm)),',',}
            round(foodavail),',',round(0.1*foodeaten),',',costofliving,',',eggslaid,',',round(eggslaid*eggmassconst));
            end;
          until (lh^[i].state='+') or (lh^[i].state='P') or (i=maxlife);
          close(invfile);
          end;

function worm.AdultAge:smallint;
         var i :smallint;
         begin
         i:=3;
         repeat inc(i) until (lh^[i].state='A') or (i=maxlife);
         if (i=maxlife) then result:=0
         else result:=age-i;
         end;
(*
function worm.AdultAge:smallint;
         var aa :smallint;
         begin
         if (age<=0) then age:=round(t-birthday-1);
         if (age<maxlife) then aa:=age else aa:=maxlife;
         repeat write(aa:6); dec(aa) until (lh^[aa].state='L') or (aa=1);
         result:=age-aa;
         end;
*)         
{$endif INDOUTPUT}

(*
procedure WriteToSampleFile;
          var i :word;
          begin
          if (samplecount>=samplen) then begin sampleon:=false; exit; end;
          inc(samplecount);
          append(samplefile);
          for i:=1 to nh do begin
            if ((h[i].modeofdeath>' ') and (sampletype='A')) or (h[i].modeofdeath='O') or ((h[i].modeofdeath>' ') and (sampletype='M') and (h[i].mature)) then begin
              WriteOneLH(i);
              inc(samplecount);
              if (samplecount>=samplen) then begin close(samplefile); sampleon:=false; exit; end;
              end;
            end;
          close(samplefile);
          end;
*)

procedure CalculateMeanAndStdev;
          var i,firstt,lastt,count          :integer;
          sum,sum2,p                        :double;
          begin
          lastt:=round(t)-2; firstt:=200; count:=succ(lastt-firstt);
          if (count<=0) then begin mean:=0; stdev:=0; exit; end;
             {compute mean and its square}
          sum:=0; sum2:=0;
          for i:=firstt to lastt do begin p:=pop[i]; sum:=sum+p; sum2:=sum2+sqr(p); end;
          mean:=sum/count;
          stdev:=sqrt(sum2/count-sqr(mean));
          end;

procedure Fourier;
          const twopi=2*pi;  maxday=20;  maxauto=200;
          var n                       :integer;
              power                   :array[1..maxday] of double;
              autocorr                :array[0..maxauto] of double;
              sqrmean                 :double;
              tempfile                :text;
              ws                      :string;

         {This isn't your standard Fourier transform. n is a number of days, and the range is truncated to be a multiple of n.
          Power is normalized by a sum over the truncated number of days}

    procedure OneSum;
              var i,firstt,lastt              :integer;
              popsum,sinsum,cossum,omega      :double;
              begin
              lastt:=round(t)-2; sinsum:=0; cossum:=0; popsum:=0;
              firstt:= 1 + (lastt mod (n*8));
              omega:=twopi/(n*8);
              for i:=firstt to lastt do begin popsum:=popsum+pop[i]; sinsum:=sinsum + sin(omega*i)*pop[i]; cossum:=cossum + cos(omega*i)*pop[i]; end;
              if (popsum=0) then power[n]:=0 {avoid divide by zero if it's a short run}
              else power[n]:=1E4*(sqr(sinsum)+sqr(cossum))/sqr(popsum);
              end;

    function  AutoCorrelation:boolean;
              var n,i,firstt,lastt,count        :integer;
              x,sum                             :double;
              begin
              result:=false;
              lastt:=round(t)-2; firstt:=200; count:=succ(lastt-firstt); if (count<1500) then exit;
              result:=true;
                {compute mean and its square}
              sum:=0;
              for i:=firstt to lastt do sum:=sum+pop[i];
              mean:=sum/count;
              sqrmean:=sqr(mean);
                 {variance is autocorr[0]}
              for n:=0 to maxauto do begin
                 sum:=0;
                 for i:=firstt to lastt-n do begin x:=pop[i]; sum:=sum+x*pop[i+n]; end;
                 for i:=lastt-n+1 to lastt do begin x:=pop[i]; sum:=sum+x*pop[i+n-lastt]; end; {wraparound}
                 autocorr[n]:=sum/count-sqrmean;
                 end;
              stdev:=sqrt(autocorr[0]+tiny);
              end;

          begin {Fourier}
          assign(fout,foutname); append(fout);
          for n:=1 to maxday do OneSum;
          writeln(fout,crlf,'Fourier analysis by number of days');
          for n:=1 to pred(maxday) do write(fout,n,','); writeln(fout,maxday);
          for n:=1 to pred(maxday) do write(fout,power[n]:5:2,','); writeln(fout,power[maxday]:5:2);
          writeln(fout,',Mean=,',mean:8:2,', StDev=,',stdev:8:2);
          if AutoCorrelation then begin
            writeln(fout,'Autocorrelation for dif=1 to ',maxauto,' normalized to variance=autocorr[0] and split into two lines:');
            for n:=1 to pred(maxauto div 2) do write(fout,autocorr[n]/autocorr[0]:6:2,','); writeln(fout,autocorr[maxauto div 2]/autocorr[0]:6:2);
            for n:=succ(maxauto div 2) to pred(maxauto) do write(fout,autocorr[n]/autocorr[0]:6:2,','); writeln(fout,autocorr[maxauto]/autocorr[0]:6:2,crlf);
            end;
          assign(tempfile,tempfilename);reset(tempfile);
          repeat
            readln(tempfile,ws); writeln(fout,ws);
          until (eof(tempfile));
          close(fout); close(tempfile);
          end; {Fourier}

procedure RateHike(which :byte; howmuch :double; howmany :integer);
          begin
          massr[which]:=massr[which]+howmuch;
          inc(countr[which],howmany);
          end;
(*
procedure UpdateAvgRates;
          var i:byte;
          begin
          for i:=1 to 15 do avcountr[i]:=lookback*avcountr[i]+oneminuslookback*countr[i];
          for i:=1 to 20 do avmassr[i]:=lookback*avmassr[i]+oneminuslookback*massr[i];
          end;
*)
{$ifdef NDEBUG}
procedure ComputeLarvamass;
          var i :integer;
          begin
          mlarvae:=0;
          for i:=1 to nworms do if (w[i].s=larva) then mlarvae:=mlarvae+w[i].leanmass {+w[i].fatmass};
          end;

procedure ComputeDauermass;
          var i :integer;
          begin
          mdauers:=0;
          for i:=1 to nworms do if (w[i].s=dauer) then mdauers:=mdauers+w[i].leanmass {+w[i].fatmass};
          end;

procedure ComputeBagmass;
          var i :integer;
          begin
          mbags:=0;
          for i:=1 to nworms do if (w[i].s=wormbag) then mbags:=mbags+w[i].leanmass {+w[i].fatmass} +w[i].eggmass;
          end;
{$endif NDEBUG}

procedure LarvaCheck;
          var i :integer;
          begin                    exit;
          mlarvae:=0; nlarvae:=0;
          if (t<100) then exit;
          for i:=1 to nworms do if (w[i].s=larva) then begin
            inc(nlarvae);
            mlarvae:=mlarvae+w[i].leanmass {+w[i].fatmass};
            end;
          if ( round(mlarvae-oldmlarvae) <> round(massr[4]+massr[6]+massr[13]-massr[7]-massr[8]-massr[9]-massr[10]-massr[11]) )
                                                         or
             ( nlarvae-oldnlarvae <> countr[4]+countr[13]-countr[8]-countr[9]-countr[10]-countr[11] )
           then begin
             write('Larva mismatch.'); readln; 
             end;
          end;

procedure ComputeStats;
          var i             :integer;
              avggen        :double;

          begin
          if (t>200) then begin
            inc(lastpd); if (lastpd<=arraysiz) then pop[lastpd]:=nworms;
            end;
          neggs:=0; nlarvae:=0; nadults:=0; ndauers:=0; nbags:=0; avggen:=0; PRcount:=0; 
                    mlarvae:=0; madults:=0; mdauers:=0; mbags:=0; {total masses in each of these phases}
          for i:=1 to nworms do with w[i] do begin
            avggen:=avggen+ gen;
            case s of
              egg     : inc(neggs);
              larva   : begin inc(nlarvae); mlarvae:=mlarvae+leanmass {+fatmass}; end;
              adult   : begin inc(nadults); madults:=madults+leanmass {+fatmass}+eggmass; if ((rage>40) and (AgeAdjustedFertility(0)<2)) then inc(PRcount); end;
              dauer   : begin inc(ndauers); mdauers:=mdauers+leanmass {+fatmass}; end;
              wormbag : begin inc(nbags); mbags:=mbags+leanmass {+fatmass}+eggmass; end;
              end;
            end;
          if (PRcount>peakPRcount) then peakPRcount:=PRcount;
          mworms:=neggs*eggmassconst+mlarvae+madults+mdauers+mbags;
          if (nworms>0) then avggen:=avggen/nworms else avggen:=0;
          LarvaCensus;
          overflow:=(nworms>maxwormmargin);
          {if (foodadded<0) then food:=-foodadded;}
          if (woutput) then begin
            {writeln(t:6:0,nworms:8, neggs:8, nlarvae:8, nadults:8, ndauers:8, nbags:8, food:8:0, starvect:8,':',RightPadItoa(maturect,8), bagct:8,':',RightPadItoa(maturect,8));}
            writeln(t:6:0,nworms:8, neggs:8, nlarvae:8, nadults:8, ndauers:8, nbags:8, mworms:12:0,eggmassconst*neggs:12:0,mlarvae:12:0,madults:12:0,mdauers:12:0,mbags:12:0, food:8:0, starvect:8,':',RightPadItoa(maturect,8));
            {$I-} append(fout); {$I+} if (ioresult>0) then rewrite(fout);
{'  time, av-gen, nworms,  neggs,nlarvae,nadults,ndauers,nbags, mworms,meggs,mlarvae,madults,mdauers,mbags,   Food,  Starve,  LarvFail,  Mature,     Bag,   DieOA,      ...Census of larvae by size... ,,,,,'+
                ' OA@D,    AA@D,    LA@D,    AM@D,    LM@D, TotPred,  Pred-E,  Pred-L,  Pred-D,  Pred-A,  Pred-B,     A@M,    T@A, EggsPerA, RepSpan,' +
                '   FoodL,   FoodA,    XrnL,    XrnA, NewEggs,  Mature,  Awaken, Brt-Dth';}
            if (upcase(foutname[length(foutname)])='V') then begin
                 write(fout,t:6:0,',',avggen:6:1,',', nworms:7,',', neggs:7,',', nlarvae:7,',', nadults:7,',', ndauers:7,',', nbags:7,',',
                 mworms:12:0,',',eggmassconst*neggs:12:0,',',mlarvae:12:0,',',madults:12:0,',',mdauers:12:0,',',mbags:12:0,',',
                 food:7:0,',', starvect:8,',',fadeawayct:10,',',maturect:8,',', bagct:8,',',oldct:8,',    ');
                 if (nlarvae=0) then write(fout,'      ,      ,      ,      ,      ,') else for i:=1 to 5 do write(fout,census[i]/nlarvae:6:3,',');
                 if (oldct=0) then oldct:=1; if (adultct=0) then adultct:=1; if (larvact=0) then larvact:=1; 
                 write(fout,olda/oldct:8:1,',',adulta/adultct:8:1,',',larvaa/larvact:8:1,',',adultm/adultct:8:1,',',larvam/larvact:8:1,',',predct:8,',',pred_e:8,',',pred_l:8,',',pred_d:8,',',pred_a:8,',',pred_b:8,',');
                 if (maturect>0) then write(fout,AatM/maturect:8:1);
                 writeln(fout,',',TasA/adultct:8:1,',',EggsPerA/AdultCt:8:1,',',RepSpan/AdultCt:8:1,',',massr[6]:8:0,',',massr[15]:8:0,',',massr[7]:8:0,',',massr[16]:8:0,',',NewEggs:8,',',MatureCt:8,',',AwakenCt:8,',',BirthsMinusDeaths:8);
                 if (round(t) mod headerepeat=headerepeat-1) then writeln(fout,csvheader);
                 end
             else
                 writeln(fout,t:6:0,nworms:8, neggs:8, nlarvae:8, nadults:8, ndauers:8, nbags:8,
                 mworms:12:0,eggmassconst*neggs:12:0,mlarvae:12:0,madults:12:0,mdauers:12:0,mbags:12:0,
                 food:8:0, (starvect+fadeawayct):8,':',RightPadItoa(maturect,8), bagct:8,':',RightPadItoa(oldct,8));
            if overflow then begin
              write('Too many worms...(press Enter)'); writeln(fout,'Too many worms'); readln;
              end;
            close(fout);
            {UpdateAvgRates;  {maintained as a moving avg, exponentially into the past}
            append(ratesfile);
            write(ratesfile,round(t):8,',');
            for i:=1 to 28 do begin
              write(ratesfile,ftoanoE(massr[i],6,1),',');
              if (i in rateswithcounts) then write(ratesfile,countr[i],',');
              end;
            {rateswithcounts:set of byte=[4,5,8,9,10,11,12,13,14,17..23];}
            {t,FI,FP,FX,EL,nEL,EP,nEP,FL,LM,LA,nLA,LD,nLD,LP,nLP,LS,nLS,DT,nDT,DL,nDL,DP,nDP,FA,AM,AA,nAA,AB,nAB,AE,nAE,AP,nAP,BB,nBB,BP,nBP,BD,nBD,'+ (BD is #23)
              'BW,LC,LZ,AC,AZ}

            {Flux identities}
           {'dF,FI-FP-FX-FL-FA,dE,AE-EL-EP,dnE,nAE-nEL-nEP,dL,EL+FL+DL-LD-LM-LP-LA-LS,dnL,EL+DL-LD-LP-LA-LS,dD,LD+BD-DL-DP-DT,dnD,LD+BD-DL-DP-DT,'+
            'dA,LA+FA-AM-AE-AP-AA-AB,dnA,LA-AP-AA-AB,dB,AB-BP-BB-BD,dnB,AB-BP-BB,Food discrep,Egg discrep,,Larva discrep,,Dauer discrep,,Adult discrep,,Bag discrep,,Waste Discrepancies';}
            write(ratesfile,round(food-oldfood),',',round(massr[1]-massr[2]-massr[3]-massr[6]-massr[15]),',');
            write(ratesfile,round((neggs-oldneggs)*eggmassconst),',',round(massr[19]-massr[4]-massr[5]),',');
            write(ratesfile,neggs-oldneggs,',',countr[19]-countr[4]-countr[5],',');
            write(ratesfile,round(mlarvae-oldmlarvae),',',round(massr[4]+massr[6]+massr[13]-massr[7]-massr[8]-massr[9]-massr[10]-massr[11]),',');
            write(ratesfile,nlarvae-oldnlarvae,',',countr[4]+countr[13]-countr[8]-countr[9]-countr[10]-countr[11],',');
            write(ratesfile,round(mdauers-oldmdauers),',',round(massr[9]+massr[23]-massr[12]-massr[13]-massr[14]),',');
            write(ratesfile,ndauers-oldndauers,',',countr[9]+countr[23]-countr[12]-countr[13]-countr[14],',');
            write(ratesfile,round(madults-oldmadults),',',round(massr[8]+massr[15]-massr[16]-massr[17]-massr[18]-massr[19]-massr[20]),',');
            write(ratesfile,nadults-oldnadults,',',countr[8]-countr[17]-countr[18]-countr[20],',');
            write(ratesfile,round(mbags-oldmbags),',',round(massr[18]-massr[21]-massr[22]-massr[23]),',');
            write(ratesfile,nbags-oldnbags,',',countr[18]-countr[21]-countr[22],','); {not 23, because that would be the count of dauers}
            {Discrepancies}
            write(ratesfile,round(food-oldfood-(massr[1]-massr[2]-massr[3]-massr[6]-massr[15])),',');
            write(ratesfile,round(eggmassconst*(neggs-oldneggs)-(massr[19]-massr[4]-massr[5])),',');
            write(ratesfile,neggs-oldneggs-(countr[19]-countr[4]-countr[5]),',');
            write(ratesfile,round(mlarvae-oldmlarvae -(massr[4]+massr[6]+massr[13]-massr[7]-massr[8]-massr[9]-massr[10]-massr[11])),',');
            write(ratesfile,nlarvae-oldnlarvae-(countr[4]+countr[13]-countr[9]-countr[8]-countr[10]-countr[11]),',');
            write(ratesfile,round(mdauers-oldmdauers -(massr[9]+massr[23]-massr[12]-massr[13]-massr[14])),',');
            write(ratesfile,ndauers-oldndauers -(countr[9]+countr[23]-countr[12]-countr[13]-countr[14]),',');
            write(ratesfile,round(madults-oldmadults -(massr[8]+massr[15]-massr[16]-massr[17]-massr[18]-massr[19]-massr[20])),',');
            write(ratesfile,nadults-oldnadults-(countr[8]-countr[17]-countr[18]-countr[20]),',');
            write(ratesfile,round(mbags-oldmbags - (massr[18]-massr[21]-massr[22]-massr[23])),',');
            write(ratesfile,nbags-oldnbags -(countr[18]-countr[21]-countr[22]),',');
            write(ratesfile,round(massr[7]-massr[25]-massr[26]),',',round(massr[16]-massr[27]-massr[28]),crlf);

            oldfood:=food;
            oldneggs:=neggs; oldnlarvae:=nlarvae; oldnadults:=nadults; oldndauers:=ndauers; oldnbags:=nbags;
                             oldmlarvae:=mlarvae; oldmadults:=madults; oldmdauers:=mdauers; oldmbags:=mbags;

            close(ratesfile);
            end; {wouotput}
{$ifdef PHEROMONE}
          UpdateCumagesum;
          inc(dauercount,ndauers);
          {adjustedfoodthreshold:=dauerfoodthresh*exp(lambda*(cumagesum/avgcumagesum - 1) );  before 8/8/16, based on cumagesum; after, based on PRcount}
          adjustedfoodthreshold:=dauerfoodthresh*exp(lambda*(PRcount/refPRcount - 1) );  {after 8/8/16, based on PRcount}
{$endif PHEROMONE}
          InitializeCountsToZero;
          end;

{$ifdef RDEBUG}
procedure Checkup;
          var i,madults,mdauers,mlarvae,mbags,meggs :integer;
          begin
          if (t<50) then exit;
          madults:=0; mdauers:=0; mlarvae:=0; mbags:=0; meggs:=0;
          for i:=1 to nworms do
            case w[i].s of adult: inc(madults);
                           dauer: inc(mdauers);
                           larva: inc(mlarvae);
                         wormbag: inc(mbags);
                             egg: inc(meggs);
                           end;
          if (frac(t+tiny)<0.0001) then exit;
          if ((madults+mdauers+mlarvae)-(nadults+ndauers+nlarvae)<>birthsminusdeaths) then begin write('Checkup BirthsMinusDeaths...'); readln; stopper:=true; end;
          if (countr[2]-countr[6]-countr[11]<>madults-nadults) then begin write('Checkup ADULTS...'); readln; stopper:=true; end;
          if (countr[1]+countr[5]-countr[4]-countr[2]-countr[9] <> mlarvae-nlarvae) then begin write('Checkup LARVAE...'); readln; stopper:=true; end;
          if (countr[4]+countr[7]-countr[5]-countr[10] <> mdauers-ndauers) then begin write('Checkup DAUERS...'); readln; stopper:=true; end;
          if (countr[6]-countr[12] <> mbags-nbags) then begin write('Checkup BAGS...'); readln; stopper:=true; end;
          if (countr[3]-countr[1]-countr[8] <> meggs-neggs) then begin write('Checkup EGGS...'); readln; stopper:=true; end;
          stopper:=false;
          end;
{$endif RDEBUG}

function worm.siz:double; {length in mm}
         begin
         result:=exp(ln(leanmass/fulladultmass)/2.7);  {mass scales with the 2.7th power of length}
         end;

function worm.stagespecificcull:double; {picks out one of the five culling input parameters}
         begin
         case s of
           egg     : result:=eggcullportion;
           larva   : result:=larvacullportion;
           dauer   : result:=dauercullportion;
           adult   : result:=adultcullportion;
           wormbag : result:=bagcullportion;
           end;
         end;

function tanh(x :double):double;
         var ex :double;
         begin
         if (x>6) then begin result:=one; exit; end;
         ex:=exp(2*x);
         result:=(ex-1)/(ex+1);
         end;

function worm.Appetite(var fertrq:double):double;
         var b,K,frac :double;
         begin
         {appetite limits lim1 and lim2 been supplanted by the logistic model for growth which is already food-dependent
  		   appetite is computed from a logistic function based on growth data (for larvae) or from growth and fertility data (for adults)}
         if (food<1000) then begin result:=small; fertrq:=0; exit; end; {prevents divide by zero and adds efficiency if there's not enough food to speak of}
         b:=bm+bn/food;  {K:=Kr+Ks/food; if (K<0) then K:=0;} K:=Kr*tanh(Ks*food);
         frac:=one - b*leanmass; if (frac<0) then frac:=0;
         result:= K*leanmass * frac / efficiency + CoL*leanmass;  {Note: I tried exp(K/8)-1 and overshot mass targets.}
         case s of
           larva : begin fertrq:=0; if (food<mworms) and (feedingdisparity>0) then result:=result*exp(feedingdisparity*ln(leanmass)); end;
           adult : begin fertrq:=AgeAdjustedFertility(4)*eggmassconst/efficiency; result:=result + fertrq; end;
           else result:=0;
           end;
         if (result<0) then result:=0;
         end;

function worm.stage:string;
         begin
         case s of
           egg : result:='egg';
           dauer : result:='dauer';
           larva : result:='larva';
           adult : result:='adult';
           wormbag : result:='wormbag';
           dead : result:='dead';
           else  result:='unknown stage';
           end
         end;

{$ifdef INDOUTPUT}
function worm.stageLetter:char;
         begin
         case s of
           egg : result:='E';
           dauer : result:='D';
           larva : result:='L';
           adult : result:='A';
           wormbag : result:='B';
           dead : result:='+';
           else  result:='?';
           end
         end;
{$endif INDOUTPUT}

procedure worm.Eat;
          var portion,share,satiety,app,p,eggmassrq :double;
          begin
{$ifdef INDOUTPUT}
          if (samplecount<samplen) and (age<=maxlife) then lh^[age].foodavail:=food;
{$endif INDOUTPUT}
          app:=Appetite(eggmassrq);
          if (eggmassrq>leanmass-minreducedmass) then
            if (leanmass>minreducedmass) then eggmassrq:=leanmass-minreducedmass else eggmassrq:=0;
          leanmass:=leanmass-eggmassrq;  eggmass:=eggmass+eggmassrq*efficiency;
          if (app<=tiny) then exit; {prevents divide by zero error}
          share:=rememberedFood*app/appsum;  {pro rata division of all the food in the tube at start of period}
          portion:= app*share /(share+app);  {limited either by pro rata share of available food, or by appetite}
          if (portion>food) then portion:=food;
          food:=food-portion;
          leanmass:=leanmass+portion*efficiency;
          satiety:=portion/app;
          totfood:=totfood+satiety;
          case s of
            larva : begin RateHike(6,portion*efficiency,0); RateHike(26,portion*(1-efficiency),0); end;
            adult : begin
                    RateHike(15,portion*efficiency,0); RateHike(28,portion*(1-efficiency),0);
                    p:=one/rage; if (p>one) then p:=one;  avgadultfood:=p*food + (one-p)*avgadultfood;
                    end;
            else begin write('You found a bug.  I shouldn''t be eating if I am not a larva or an adult.  Press <Enter>...'); readln; end;
            end;
{$ifdef INDOUTPUT}
          if (samplecount<samplen) and (age<=maxlife) then lh^[age].foodeaten:=round(portion*10);
{$endif INDOUTPUT}
          end;

function worm.rage:double;
         const egg2dauertime=8;
         begin
         result:=t-ttime+small;
         if (s=larva) and (enterdauer_t>0) then result:=result + (enterdauer_t-birthday); {rage for a larva is age excluding time spent as dauer}
         end;

function worm.AgeAdjustedFertility(setback :double):double;  {"setback" is 4/8 of a day because measured fertilities are for trailing day.  }
         var aage,ri,rj :double;
             i,j        :byte;

   function Interpolate:double;
            var p :double;
            begin
            p:=(avgadultfood-fd[i]) / (fd[j]-fd[i]);
            result:=ri*(1-p)+rj*p;
            if (result<0) then result:=0;
            end;

   function InterpolateLogarithmically:double;
            var p :double;
            begin
            p:=(avgadultfood-fd[i]) / (fd[j]-fd[i]);
            result:=exp(ln(ri)*(1-p)+ln(rj)*p);
            end;

         begin
         aage:=(rage+setback)/8/fertspan;
         if (avgadultfood>=fd[5]) then begin result:=fertmult[5]*exp(xp[5]*ln(aage)-aage/agescale[5]); exit; end;
         if (avgadultfood<fd[2]) then i:=1 else if (avgadultfood<fd[3]) then i:=2 else if (avgadultfood<fd[4]) then i:=3 else i:=4;
         j:=i+1;
         ri:=fertmult[i]*exp(xp[i]*ln(aage)-aage/agescale[i]);
         rj:=fertmult[j]*exp(xp[j]*ln(aage)-aage/agescale[j]);
         result:= InterpolateLogarithmically {* leanmass/fulladultmass};
         end;

procedure TestFertility;
          var w    :worm;
              i,j  :word;
              f    :text;
          begin
          assign(f,'FertTest.csv');rewrite(f);
          w.leanmass:=fulladultmass;
          t:=200;
          for j:=1 to 6 do begin
            w.avgadultfood:=fd[j];
            for i:=0 to 160 do begin
              w.ttime:=t-i;
              write(f,w.AgeAdjustedFertility(4):6:2,',');
              end;
            writeln(f);
            end;
          close(f);
          end;

procedure worm.Hatch;
          begin
          if (synchronized) and (t>25) then exit;
          if (rage>hatchage) or ((4*rage>3*hatchage) and (random<quarter)) then begin  {independent of everything.  hatchage is now fixed, not user-specified}
            inc(allworms);
            RateHike(4,eggmassconst,1);
            inc(BirthsMinusDeaths);
            s:=larva;
            ttime:=t;
            age:=0;
            totfood:=half*foodunit; {don't initialize this to zero or else the worm will never eat}
            {fatmass:=half*newlarvamass;
            leanmass:=half*newlarvamass;}
            if (food>newlarvamass) then begin
  			      leanmass:=newlarvamass;
              food:=food - (newlarvamass-eggmassconst);
              RateHike(6,newlarvamass-eggmassconst,0);
              end;
            eggmass:=0;
            lastfood:=food;
{$ifdef VARIABILITY}
            ipumplimit:=exp(GaussRand*indvariation)*pumplimit;
            iminmaturemass:=exp(GaussRand*indvariation)*maturemassmin;
            imaturityage:=exp(GaussRand*indvariation)*minmaturityage;
{$endif VARIABILITY}
            end;
          {$ifdef RDEBUG} Checkup; {$endif RDEBUG}
          end;

procedure worm.Mature;  {Only get here if age is qualified for maturity} 
         var lmass,amass,
             exrn,delta   :double;
          begin
          if (leanmass<minmaturemass) and (rage>maxmaturityage) then begin Starve; exit; end;
          if (synchronized) and (t>20+minmaturityage) then begin Starve; exit; end;
          if (leanmass>minmaturemass) then begin
            inc(matureworms);
             {Old formula derived algebraically from three equations: Let delta equal mass lost to anabolic inefficiency.  then
                       fatmass1 + leanmass1 + delta = fatmass0 + leanmass0
                       fatmass1 = leanmass1 * adult_app_c0
                       delta = (1-efficiency) * (fatmass0-fatmass1)                         }
             {6/17 formula derived algebraically from two equations:
                       fatmass1=fatatmaturity*leanmass1   and
                       leanmass1=[leanmass0 + (fatmass0-fatmass1)] * efficiency             }
             {9/17 - Don't reduce lean mass; fatamaturity=1/6 is a maximum only, and fat can be less than 1/6}
            minreducedmass:= leanmass*(1-leanmass2eggsproportion) + minmaturemass*leanmass2eggsproportion; 
            lmass:=leanmass;
            avgadultfood:=food;
            exrn:=lmass*CoL;
            delta:=leanmass-lmass-exrn; {efficiency loss}
            amass:=leanmass;
            RateHike(7,exrn+delta,0);
            RateHike(25,exrn,0);
            RateHike(26,delta,0);
            RateHike(8,leanmass,1);
            inc(maturect);
            AatM:=AatM+rage;
            s:=adult;
            ttime:=t;
{$ifdef INDOUTPUT}
            if (samplecount<samplen) then matur:=true;
{$endif INDOUTPUT}
            end;
          {$ifdef RDEBUG} Checkup; {$endif RDEBUG}
          end;

procedure worm.Alchemy;
          var x,exrn,fat2lean,waste :double;
          begin
          if (s<>larva) and (s<>adult) then exit;
          exrn:=leanmass*CoL;
          if (s=larva) then begin
             RateHike(7,exrn,0);
             RateHike(25,exrn,0);
             end
           else begin
             RateHike(16,exrn,0);
             RateHike(27,exrn,0);
             end;
          leanmass:=leanmass - exrn;  {cost of living deduction}
          end; {Alchemy}

function worm.Attrition:boolean; {Even dauers die after a time}
         const m=6.0;
         var ex,aa :double;
         begin
         aa:=rage;
         ex:=exp(m*(rage/attrition_time - 1));
         result:=(random<ex/(1+ex));
         if result then begin
           RateHike(12,leanmass{+fatmass},1);
           s:=dead;
           dec(BirthsMinusDeaths);
           TasA:=TasA+rage;
{$ifdef INDOUTPUT}
           if (samplecount<samplen) then begin
             if (age<=maxlife) then lh^[age].state:='+';
             modeofdeath:='T';
             end;
{$endif INDOUTPUT}
           end;
         end;

procedure worm.Awaken;
          begin
{$ifdef INDOUTPUT}
          if (samplecount<samplen) and (age<=maxlife) then lh^[age].foodavail:=food;
{$endif INDOUTPUT}
          {if (food<dauerfoodthresh) or (random>dauerawakenprob) then exit;  THIS IS THE CONDITIONAL OPTION, NOT AS GOOD 2/12/16}
{$ifdef PHEROMONE}
          {if (random<dauerawakenprob * (1-exp(-(food+lastfood)/(adjustedfoodthreshold)))) then begin Replaced by sqrt model from lab data, 12018Dec}
		  if (random<dauerawakenprob * sqrthalf * sqrt(half*(food+lastfood)*dauerfoodthreshold/adjustedfoodthreshold) then begin
{$else PHEROMONE}
          {if (random<dauerawakenprob * (1-exp(-(food+lastfood)/(dauerfoodthresh)))) then begin  Replaced by sqrt model from lab data, 12018Dec}
		      if (random<dauerawakenprob * sqrt(half*(food+lastfood))) then begin
{$endif PHEROMONE}
            s:=larva; ttime:=t; 
            inc(awakenct);
            RateHike(13,leanmass{+fatmass},1);
{$ifdef INDOUTPUT}
          if (samplecount<samplen) then begin
              if (age<=maxlife) then lh^[age].state:='W';
              end;
{$endif INDOUTPUT}
            end;
          {$ifdef RDEBUG} Checkup; {$endif RDEBUG}
          end;

procedure worm.NewWorm(gg :integer);
          var j :word;
          begin
          age:=0;
          gen:=gg;
          s:=egg;
          cum_eggs:=0;
          frombag:=false; {will be changed if it's a bagger}
          birthday:=t*(1+random*tiny); {tiny random is to make sure birthdays are unique}
          tatfirstegg:=0; tatlastegg:=0;
          ttime:=birthday;
          enterdauer_t:=0;
          totfood:=0;
          lastfood:=food;
          leanmass:=dauermassconst; {fatmass:=0;} eggmass:=0;  {these last three will be overwritten}
{$ifdef INDOUTPUT}
            if (lh=nil) then lh:=new(tsptr); 
            modeofdeath:=' '; matur:=false;
            for j:=1 to maxlife do with lh^[j] do begin
              {leanm:=999; fatm:=999; foodeaten:=999; eggslaid:=99;}
              leanm:=0; {fatm:=0;} foodeaten:=0; eggslaid:=0;
              state:='X';
              end;
            with lh^[0] do begin
              state:='E';
              leanm:=round(half*eggmassconst);
              {fatm:=leanm;}
              foodavail:=food;
              foodeaten:=0;
              eggslaid:=0;
              end;
{$endif INDOUTPUT}
          end;

procedure worm.Reproduce;
{Every egg gets laid when it is formed.}
          var i,fert :word;
          begin
          if (eggmass<eggmassconst) then exit;
(*
          if (eggmass/eggmassconst>100) then begin
            write(round(eggmass/eggmassconst),' Dangerous call to Poisson');
            fert:=round(eggmass/eggmassconst * (1 + 0.1*(random-random)));
            end                                                                                                               f
          else fert:=FastPoissonDist(eggmass/eggmassconst);
*)
          fert:=floor(eggmass/eggmassconst);  {Kerry's request, 6/13/16}
          if (fert=0) then exit;
          tatlastegg:=t;
          if (cum_eggs=0) then tatfirstegg:=t;
          inc(cum_eggs,fert);
          inc(neweggs,fert);
          eggmass:=eggmass-fert*eggmassconst;
          RateHike(19,fert*eggmassconst,fert);
          for i:=1 to fert do if (nworms<maxwormmargin) then begin
            inc(nworms);
            w[nworms].NewWorm(succ(gen));
            end;
{$ifdef INDOUTPUT}
          if (samplecount<samplen) then begin
            if (age<=maxlife) then lh^[age].eggslaid:=fert;
            end;
{$endif INDOUTPUT}
          {$ifdef RDEBUG} Checkup; {$endif RDEBUG}
          end;

function worm.Bag:boolean;
         begin
         result:=false;
         {if (s<>adult) or (fatmass>0) or (random>AgeAdjustedFertility) then exit;}
         if (s<>adult) then exit;
(*         if (fatmass<tiny) then begin
         HERE's THE PROBLEM (2018Nov)
            if ((rage>40) and (AgeAdjustedFertility(0)<2)) then begin Die; exit; end; {lumped with dying of old age}
            end; *)
         if (food+lastfood>bagthreshold) then exit;
         result:=true;
         TasA:=TasA+rage;
         {to be computed on bagging eggcount:=trunc((leanmass+fatmass+eggmass)*larva_anaboleffic*twothirds/dauermassconst); {New criterion for bagn 2016Aug18}
         {eggmass:=eggmassconst*neggs; leave eggmass alone. For bags it will not agree with the number of eggs, but the total mass of bags will be accurate.}
         EggsPerA:=EggsPerA+cum_eggs+neggs;
         inc(NewEggs,neggs);
         RepSpan:=RepSpan+tatlastegg-tatfirstegg;
         s:=wormbag;
         ttime:=t;
         inc(bagct);
         inc(adultct);
         dec(BirthsMinusDeaths);
         adulta:=adulta+rage;
         adultm:=adultm+leanmass+eggmass{+fatmass};
         RateHike(18,leanmass{+fatmass}+eggmass,1);
{$ifdef INDOUTPUT}
         if (samplecount<samplen) then begin
           modeofdeath:='B';
           if (age<=maxlife) then lh^[age].state:='B';
           end;
{$endif INDOUTPUT}
         {$ifdef RDEBUG} Checkup; {$endif RDEBUG}
         end;

procedure worm.Burst;
          var i,bagn :word;
              waste  :double;
          begin
          if (t<ttime+burstage) then exit;
          {bagn:=round(eggmass/eggmassconst);  {eggmass, and thus bagn was determined at the time of bagging}
          {bagn:=trunc((leanmass+fatmass+eggmass)*efficiency*twothirds/dauermassconst); {New criterion for bagn 2016Aug18}
          {if (bagn>neggs) then bagn:=neggs;}
          bagn:=trunc((leanmass{+fatmass}+eggmass)*efficiency*twothirds/dauermassconst); {New criterion for bagn 2016Aug18}
          waste:=leanmass{+fatmass}+eggmass-bagn*dauermassconst;  
          RateHike(23,bagn*dauermassconst,bagn);   {create the dauers}
          RateHike(21,waste,1);  {destroy the bag}
          {RateHike(24,waste,1); deleted because all accounted for in 21, as of 2017Oct}
          for i:=1 to bagn do if (nworms<maxwormmargin) then begin
            inc(nworms);
             {Don't use "with" because of ambiguity when placed within an object that refers to another instance of worm}
            w[nworms].NewWorm(succ(gen));
            w[nworms].frombag:=true;
            w[nworms].age:=0;
            w[nworms].s:=dauer;
            w[nworms].birthday:=t-(burstage div 2)+tiny*(random-random); {2017Aug birthday arbitrarily set at 1/2 burstage time units before bag bursts}
            w[nworms].totfood:=half*foodunit;
            w[nworms].leanmass:=dauermassconst;
            {w[nworms].fatmass:=half*dauermassconst;}
            w[nworms].lastfood:=food;
            w[nworms].enterdauer_t:=t;
{$ifdef INDOUTPUT}
            with w[nworms].lh^[0] do begin
              state:='D'; {indicating a birth from bag}
              leanm:=round(w[nworms].leanmass);
              {fatm:=round(w[nworms].fatmass);}
              foodeaten:=0;
              eggslaid:=0;
              end;
{$endif INDOUTPUT}
            end;
          inc(BirthsMinusDeaths,bagn);
          s:=dead; leanmass:=0; {fatmass:=0;} eggmass:=0;
          {$ifdef RDEBUG} Checkup; {$endif RDEBUG}
          end;

{$ifdef PHEROMONE}
function worm.Gompertz:double;
         const ef=8; {number of Gompertz e-foldings in one lifespan}
         var ls,p :double;
         begin
         {ls:=lifespan*exp(mu*(cumagesum/avgcumagesum - 1));}
{$ifdef VARIABILITY}
         tau:=tummult*ilifespan/GompertzN;
         A:=ilifespan*(exp(GompertzN)-one}
{$endif VARIABILITY}
         {result:=exp(-exp(ef*((rage+dietmult*totfood)/ls - 1)));  SUPPLANTED 1/19}
         p:=(exp(rage/tau)-one)/A; {Probability of death per unit time}
         result:=exp(-p);  {Probability of surviving one time period}
         end;
{$else PHEROMONE}
function worm.Gompertz:double;
         const ef=8; {number of Gompertz e-foldings in one lifespan}
         var p :double;
         begin
{$ifdef VARIABILITY}
         tau:=tummult*ilifespan/GompertzN;
          A:=ilifespan*(exp(GompertzN)-one}
{$endif VARIABILITY}
          p:=(exp(rage/tau)-one)/A; {Probability of death per unit time}
          result:=exp(-p);  {Probability of surviving one time period}
         end;
{$endif PHEROMONE}

procedure worm.SenescentMortality;
          var a,g :double;
          begin
          {a:=rage; g:=Gompertz;}
          if (random>Gompertz) then
               Die;
          end;

procedure worm.Die; {of old age}
          var i :word;
             {p=exp(ef*(age/lifespan-1)) is the prob per unit time of dying. 1-exp(-pdt) is probability of surviving this time step.  If time is measured in units of a time step dt, then
              the probability of death is 1-exp(-exp(age/ls)) and that is implemented as random>exp(-exp(age/ls)) }
          begin
          inc(oldageworms);
          s:=dead;
          EggsPerA:=EggsPerA+cum_eggs;
          RateHike(17,leanmass{+fatmass}+eggmass,1);
          inc(oldct); inc(adultct); dec(BirthsMinusDeaths);
          TasA:=TasA+rage;
          olda:=olda+rage;
          adultm:=adultm+leanmass{+fatmass}+eggmass;
          RepSpan:=RepSpan+tatlastegg-tatfirstegg;
{$ifdef INDOUTPUT}
          if (samplecount<samplen) then begin
            if (age<=maxlife) then lh^[age].state:='+';
            modeofdeath:='A';
            end;
{$endif INDOUTPUT}
          {$ifdef RDEBUG} Checkup; {$endif RDEBUG}
          end;

procedure worm.Starve;  {only for larvae; adults bag instead}
          begin
          if (random<escapestarvationprob) then begin
             leanmass:=leanmass*0.9;
             exit;
             end;
          inc(starvect);
          s:=dead;
          RateHike(11,leanmass{+fatmass},1);
          inc(larvact); dec(BirthsMinusDeaths);
          larvam:=larvam+leanmass{+fatmass};
          larvaa:=larvaa+rage;
{$ifdef INDOUTPUT}
          if (samplecount<samplen) then begin
            if (age<=maxlife) then lh^[age].state:='+';
            modeofdeath:='S';
            end;
{$endif INDOUTPUT}
          {$ifdef RDEBUG} Checkup; {$endif RDEBUG}
          end;

function worm.Punt:boolean;  {larva becomes a dauer}
          const tol=1.25;
          begin
             {OLD: if (leanmass/dauermassconst<tol) and (dauermassconst/leanmass<tol) and (fatmass/leanmass<larvadauerthresh) then begin s:=dauer;  ttime:=t; end;}
          {Larvae go dauer with a probability=larvadauerthresh based on 2 successive readings of food concentration during the size window}
          result:=false;
          {if (leanmass>0.6*dauermassconst) and (leanmass<2*dauermassconst) and (random<maxdaueringrate*exp(-half*(food+lastfood)/dauerfoodthresh)) then begin}
          if (random<maxdaueringrate*exp(-half*(food+lastfood)/dauerfoodthresh)) then
            if (leanmass<0.6*dauermassconst) or (leanmass>2*dauermassconst) or (enterdauer_t<>0) {mass out of dauering range or previous dauer experience} then begin result:=true; Starve; end  {Result=true is a kluge to keep this worm from eating when it's dead}
            else begin
            s:=dauer;  ttime:=t;
            enterdauer_t:=t;
            RateHike(9,leanmass{+fatmass},1);
            result:=true;
{$ifdef INDOUTPUT}
            if (samplecount<samplen) and (age<=maxlife) then with lh^[age] do begin
               state:='D';
               {if (foodadded<0) then food:=-foodadded;}
               foodavail:=food;
               leanm:=round(leanmass); {fatm:=round(fatmass);}
               end;
{$endif INDOUTPUT}
            end
          else lastfood:=food;
          {$ifdef RDEBUG} Checkup; {$endif RDEBUG}
          end;

procedure worm.Eliminate;
          var temp :worm;
          begin
{$ifdef INDOUTPUT}
          if (sampleon) and (samplecount<samplen) then
            if (sampletype in ['A','E']{All or Early}) or (modeofdeath='A' {Aging}) or ((sampletype='F') and (tatlastegg-tatfirstegg>12))
            or ((modeofdeath>' ') and (sampletype='M') and (matur)) then begin
              append(samplefile);
              WriteOneLH; {includes both individual file and line of samplefile}
              close(samplefile);
              inc(samplecount); if (samplecount>=samplen) then begin sampleon:=false; end;
              end;
          {dispose(lh); lh:=nil;}
{$endif INDOUTPUT}
          temp:=w[nworms];
          w[nworms]:=self;  {So we don't have to be allocating memory on the fly; this keeps the lh in circulation}
          self:=temp;
          dec(nworms);
          end;

procedure worm.Predate;
          begin
          inc(predct);
          case s of
              larva   : begin
                        inc(larvact); dec(BirthsMinusDeaths);
                        larvaa:=larvaa+rage; larvam:=larvam+leanmass{+fatmass};
                        inc(pred_l);
                        RateHike(10,leanmass{+fatmass},1);
                        end;
              adult   : begin
                        inc(adultct); dec(BirthsMinusDeaths);
                        RepSpan:=RepSpan+tatlastegg-tatfirstegg;
                        EggsPerA:=EggsPerA+cum_eggs;
                        TasA:=TasA+rage;
                        adulta:=adulta+rage; adultm:=adultm+leanmass{+fatmass}+eggmass;
                        inc(pred_a);
                        RateHike(20,leanmass{+fatmass}+eggmass,1);
                        end;
              wormbag : begin inc(pred_b); RateHike(22,leanmass{+fatmass}+eggmass,1); end;
              egg     : begin inc(pred_e); RateHike(5,eggmassconst,1); end;
              dauer   : begin inc(pred_d); RateHike(14,leanmass{+fatmass},1); dec(BirthsMinusDeaths); end;
              dead    : begin write('You found a bug. Dead worm where there shouldn''t be one.  Press <Enter> to continue...'); readln; end;
              end;
{$ifdef INDOUTPUT}
           if (samplecount<samplen) then begin
            if (modeofdeath=' ') then modeofdeath:='P'; {if it hasn't already bagged}
            if (age<=maxlife) then lh^[age].state:='P';
            end;
{$endif INDOUTPUT}
          s:=dead;
          Eliminate;
          {$ifdef RDEBUG} Checkup; {$endif RDEBUG}
          end;

{$ifdef INDOUTPUT}
procedure worm.RecordLHState;
          begin
          if (age>maxlife) then exit; {avoid range check errors}
          lh^[age].leanm:=round(leanmass);
          {if (fatmass<0) then lh^[age].fatm:=0
          else if (fatmass>high(word)) then lh^[age].fatm:=high(word)
          else lh^[age].fatm:=round(fatmass);}
          lh^[age].state:=stageLetter;
          end;
{$endif INDOUTPUT}

procedure FeedAndCull;
          var i :integer;
          begin
          i:=1;
          while (i<=nworms) do
            if (random<cullportion) then w[i].Predate
            else if (random<w[i].stagespecificcull) then w[i].Predate else inc(i);
          RateHike(2,food*cullportion,0);
          food:=food*(1-cullportion)*(1-foodcullportion);
          food:=food+foodadded; lastfeeding:=round(t); lastculling:=round(t); {the "round" makes sure that lastfeeding doesn't slip from a cycle based on integral number of days in the case of 7322 Mon-Wed-Fri schedule}
          RateHike(1,foodadded,0);
          if (randomfeeding) then feedinginterval:=maxterval+random(succ(maxterval-minterval));
          end;

procedure Cull;
          var i :integer;
          begin
          i:=1;
          while (i<=nworms) do
            if (random<cullportion) then w[i].Predate
            else if (random<w[i].stagespecificcull) then w[i].Predate else inc(i);
          food:=food*(1-cullportion)*(1-foodcullportion);
          RateHike(2,food*cullportion,0);
          lastculling:=round(t);
          if (randomculling) then cullinginterval:=maxcerval+random(succ(maxcerval-mincerval));
          end;

procedure Feed;
          begin
          food:=food+foodadded; lastfeeding:=round(t); {the "round" makes sure that lastfeeding doesn't slip from a cycle based on integral number of days in the case of 7322 Mon-Wed-Fri schedule}
          RateHike(1,foodadded,0);
          lastfeeding:=round(t);
          if (randomfeeding) then feedinginterval:=maxterval+random(succ(maxterval-minterval));
          end;

procedure DeListDeadWorms;
          var i :integer;
          begin
          i:=1;
          while (i<=nworms) and (i>0) do
            if (w[i].s=dead) then w[i].Eliminate else inc(i);
          end;

procedure OnePeriod;
          const changefood:integer=20000000;
          var i,j,nworms_      :integer;
              ta,interval      :double;

   procedure SumAllAppetites;
             var i  :integer;
                 ta :double;
             begin
             if (foodadded<0) then food:=-foodadded;
             appsum:=0;
             for i:=1 to nworms do if (w[i].s in [larva,adult]) then appsum:=appsum+w[i].appetite(ta);
             rememberedFood:=food;
             end;

(* { ifdef INDOUTPUT}
    procedure UpdateMasses;
              var j,age :integer;
                  i     :word;
              begin
              for j:=1 to nworms do with w[j] do begin
                {write(j:4);}
                if (age>maxlife) then continue; {avoid range check errors}
                {if (i>nh) or (i<1) then begin write('UpdateMasses: bad worm index=',i); readln; end;
                if (j>nworms) or (j<1) then begin write('UpdateMasses: bad worm number=',j); readln; end;
                write(',',i:4); write(age:4);}
                lh[age].leanm:=round(leanmass);
                if (fatmass<0) then lh[age].fatmass:=0
                else lh[age].fatm:=round(fatmass);
                end;
              end;
{ endif INDOUTPUT}
*)
   procedure Shuffle;
             var i,m     :integer;
                 swap    :worm;
             begin
             for i:=1 to nworms_ do begin
                m:=succ(random(nworms_));
                swap:=w[i];
                w[i]:=w[m];
                w[m]:=swap;
                end;
             end;

          begin {OnePeriod}
(*
if (changefood<>0) and (t>22) then begin
   foodadded:=-abs(changefood);
   changefood:=0;
   end;
*)
          if (feedinginterval=-1) then begin
            if (round(t)>=nextfeeding) then begin
              if (feedandculltogether) then FeedAndCull else Feed;
              ReadFromSchedFile;
              end;
            end
          else if (feedandculltogether) then begin
            if (feedinginterval=7223) or (feedinginterval=7322) then
              if ( (round(lastfeeding) mod 56 = 0) and (round(t-lastfeeding)>=24) )
                                          or
                 ( (round(lastfeeding) mod 56 > 0) and (round(t-lastfeeding)>=16) )
                 then FeedAndCull;
            if (round(t-lastfeeding)>=feedinginterval) then FeedAndCull; {if feedinginterval=7223 we'll never get here}
            end
          else begin
            if (cullinginterval=7223) or (cullinginterval=7322) then
              if ( (round(lastculling) mod 56 = 0) and (round(t-lastculling)>=24) )
                                          or
                 ( (round(lastculling) mod 56 > 0) and (round(t-lastculling)>=16) )
                 then Cull;
            if (round(t-lastculling)>=cullinginterval) then Cull; {if cullinginterval=7223 we'll never get here}
            if (feedinginterval=7223) or (feedinginterval=7322) then
              if ( (round(lastfeeding) mod 56 = 0) and (round(t-lastfeeding)>=24) )
                                          or
                 ( (round(lastfeeding) mod 56 > 0) and (round(t-lastfeeding)>=16) )
                 then Feed;
            if (round(t-lastfeeding)>=feedinginterval) then Feed; {if feedinginterval=7223 we'll never get here}
            end;
          if (nworms=0) then exit;
          nworms_:=nworms;
          interval:=one/nworms;
          SumAllAppetites;
          if (round(t)=100) and (sampletype<>'E') then samplecount:=0 {initializes record-keeping} else if (round(t)=200) then sampleon:=true; {initializes record-reporting}
          Shuffle;
          for i:=1 to nworms_ do begin
            t:=t+interval;
            j:=i;
            inc(w[j].age);
            if (foodadded<0) then food:=-foodadded;
            with w[j] do case s of  {each of these procedures has probabilities built in}
               egg    : Hatch;
               larva  : begin if Punt then continue; Eat; Alchemy; if (rage>=minmaturityage) then Mature; {includes tests for Starve;} end;
               dauer  : if not Attrition then Awaken;
               adult  : begin Eat; Alchemy; Reproduce; if (not Bag) then SenescentMortality; end;
               wormbag: Burst;
               dead   : begin writeln('Dead worm found where there should not be one.'); readln; end;
               end;
            {$ifdef RDEBUG} if (t>100) and (t<110) then begin writeln(t:6:3,j:8); Checkup; end; {$endif RDEBUG}
            w[j].lastfood:=food;
            appsum:=appsum-w[j].appetite(ta); if (appsum<one) then appsum:=one;
{$ifdef INDOUTPUT}
            with w[j] do if (samplecount<samplen) and (age<=maxlife) then RecordLHState;
{$endif INDOUTPUT}
            end;
          DeListDeadWorms; {invokes Eliminate, where sample output is written to disk}
          RateHike(3,food*(1-foodattrition),0);
          food:=foodattrition*food;  {food loss, according to an input variable}
          {ComputeLarvaMass;
          if (mlarvae-oldmlarvae -(massr[4]+massr[13]-massr[9]+massr[6]-massr[7]-massr[8]-massr[10]-massr[11])>1E-5) then begin write(' End of Time Step'); readln; end;}
          ComputeStats;
          end;  {OnePeriod}

procedure MainLoop;
          var i1,i2,i3,i4 :byte;

   procedure ProcessFoutName;
             var p :byte;
             begin
             foutname:=barefoutname;
             if (LL[1].max=1) then exit;
             p:=pos('.',foutname);
             {Insert('-'+itoa(i1),foutname,p);}
             Insert(itoa(i1),foutname,p);
             if (LL[2].max=1) then exit;
             p:=pos('.',foutname);
             Insert(itoa(i2),foutname,p);
             if (LL[3].max=1) then exit;
             p:=pos('.',foutname);
             Insert(itoa(i3),foutname,p);
             if (LL[4].max=1) then exit;
             p:=pos('.',foutname);
             Insert(itoa(i4),foutname,p);
             end;

          begin {MainLoop}
          for i1:=1 to LL[1].max do begin
            with LL[1] do if (max>1) then if (isdouble) then dbladr^:=valu[i1] else intadr^:=round(valu[i1]);
            for i2:=1 to LL[2].max do begin
              with LL[2] do if (max>1) then if (isdouble) then dbladr^:=valu[i2] else intadr^:=round(valu[i2]);
              for i3:=1 to LL[3].max do begin
                with LL[3] do if (max>1) then if (isdouble) then dbladr^:=valu[i3] else intadr^:=round(valu[i3]);
                for i4:=1 to LL[4].max do begin
                  with LL[4] do if (max>1) then if (isdouble) then dbladr^:=valu[i4] else intadr^:=round(valu[i4]);
                  Init;
                  ProcessFoutName;
                  OpenOutputFile(i1,i2,i3,i4);
                  try
                    repeat OnePeriod until (nworms<5) or (t+tiny>=stoptime) or (overflow) or (done);
                    Fourier;
                  except
                    on EZeroDivide do continue;
                    on EOverflow do continue;
                    on EMathError do continue;
                    end;
                  end;
                end;
              end;
            end;
          end; {MainLoop}

{$ifdef PHEROMONE}
procedure MuLoop;
          var d,imu :integer;
          begin
          ReadInputFile;
          mu:=0;
          lambda:=2; {optimum value from LambdaLoop}
          assign(pherofile,'PheroSum.csv'); {$I-}append(pherofile);{$I+} if (ioresult>0) then rewrite(pherofile);
          writeln(pherofile);

          writeln(pherofile,DateTimeToStr(date+time),'  Life span modulated by PR worms.');
          writeln(pherofile,crlf,'end-t, life span, mu, total-dauers, end-dauers, mean, stdev');

          for d:=0 to 10 do begin
            writeln(d,lifespan:8);
            for imu:=-10 to 150 do begin
              mu:=imu/50;
              adjustedfoodthreshold:=dauerfoodthresh;
              Init;
              if (woutput) then OpenOutputFile(1,1,1,1);
              repeat OnePeriod until (ndauers=nworms) or (t+tiny>=stoptime) or (overflow);
              CalculateMeanAndStdev;
              writeln(pherofile,round(t),',',lifespan,',',mu:8:3,',',dauercount,',',ndauers,',',mean:8:1,',',stdev:8:2);
              end;
            lifespan:=lifespan+5;
            end;
          close(pherofile);
          end;

procedure OneLambdaRun(savelambda,dawake :double);
          begin
          mu:=0; dauerawakenprob:=dawake;
          lambda:=0; refPRcount:=maxint;
          assign(pherofile,'PheroSum.csv'); {$I-}append(pherofile);{$I+} if (ioresult>0) then rewrite(pherofile);
          writeln(pherofile);

          writeln(pherofile,DateTimeToStr(date+time),'  Enter and exit dauer phase modulated by recent number of senior worms.');
          writeln(pherofile,crlf,'end-t, d-awake-prob, lambda, total-dauers, end-dauers');

          adjustedfoodthreshold:=dauerfoodthresh;
          Init;
          if (woutput) then OpenOutputFile(1,1,1,1);
          repeat OnePeriod until (ndauers=nworms) or (t+tiny>=stoptime) or (overflow);
          writeln(pherofile,round(t),',',dauerawakenprob:8:3,',',lambda:8:3,',',dauercount:8,',',ndauers:8);
          lambda:=savelambda; refPRcount:=round(exp(-1)*peakPRcount);
          Init;
          if (woutput) then OpenOutputFile(1,1,1,1);
          repeat OnePeriod until (ndauers=nworms) or (t+tiny>=stoptime) or (overflow);
          writeln(pherofile,round(t),',',dauerawakenprob:8:3,',',lambda:8:3,',',dauercount:8,',',ndauers:8);
          close(pherofile);
          end;

procedure OneMuRun(savemu,dawake :double);
          begin
          mu:=0; dauerawakenprob:=dawake;
          lambda:=0; refPRcount:=maxint;
          assign(pherofile,'PheroSum.csv'); {$I-}append(pherofile);{$I+} if (ioresult>0) then rewrite(pherofile);
          writeln(pherofile);

          writeln(pherofile,DateTimeToStr(date+time),'  Enter and exit dauer phase modulated by recent number of senior worms.');
          writeln(pherofile,crlf,'end-t, d-awake-prob, mu, total-dauers, end-dauers');

          adjustedfoodthreshold:=dauerfoodthresh;
          Init;
          if (woutput) then OpenOutputFile(1,1,1,1);
          repeat OnePeriod until (ndauers=nworms) or (t+tiny>=stoptime) or (overflow);
          writeln(pherofile,round(t),',',dauerawakenprob:8:3,',',mu:8:3,',',dauercount:8,',',ndauers:8);
          mu:=savemu; refPRcount:=round(exp(-1)*peakPRcount);
          Init;
          if (woutput) then OpenOutputFile(1,1,1,1);
          repeat OnePeriod until (ndauers=nworms) or (t+tiny>=stoptime) or (overflow);
          writeln(pherofile,round(t),',',dauerawakenprob:8:3,',',mu:8:3,',',dauercount:8,',',ndauers:8);
          close(pherofile);
          end;

procedure LambdaLoop;
          var d,ilambda :integer;
          begin
          ReadInputFile;
          lambda:=0;
          assign(pherofile,'PheroSum.csv'); {$I-}append(pherofile);{$I+} if (ioresult>0) then rewrite(pherofile);
          writeln(pherofile);

          writeln(pherofile,DateTimeToStr(date+time),'  Enter and exit dauer phase modulated by recent number of senior worms.');
          writeln(pherofile,crlf,'end-t, d-awake-prob, lambda, rnd(lambda), total-dauers, end-dauers, mean, stdev');
          dauerawakenprob:=0.02;
          for d:=1 to 3 do begin
            writeln(d,dauerawakenprob:8:3);
            for ilambda:=0 to 250 do begin
              lambda:=ilambda/50;
              adjustedfoodthreshold:=dauerfoodthresh;
              Init;
              if (woutput) then OpenOutputFile(1,1,1,1);
              repeat OnePeriod until (ndauers=nworms) or (t+tiny>=stoptime) or (overflow);
              CalculateMeanAndStdev;
              writeln(pherofile,round(t),',',dauerawakenprob:8:3,',',lambda:8:3,',',0.1*round(lambda*10):5:1,',',dauercount,',',ndauers,',',mean:8:1,',',stdev:8:2);
              end;
            dauerawakenprob:=dauerawakenprob*2;
            end;
          close(pherofile);
          end;

procedure PheroSurvivalLoop;
          var k,c,d,ilambda,sumt :integer;
          begin
          feedfilename:='Ever_Sparser0.inp';
          ReadInputFile;
          lambda:=0;
          assign(pherofile,'SurviveTime.csv'); {$I-}append(pherofile);{$I+} if (ioresult>0) then rewrite(pherofile);
          writeln(pherofile);

          writeln(pherofile,DateTimeToStr(date+time),'  Absolute expiry date for dauers added.  Enter and exit dauer phase modulated by recent number of senior worms.');
          writeln(pherofile,crlf,'lambda, end-t, d-awake-prob, rnd(lambda), dauercount, mean, stdev');
          dauerawakenprob:=0.10;
(*
          adjustedfoodthreshold:=dauerfoodthresh;
          Init;
          if (woutput) then OpenOutputFile(1,1,1,1);
          repeat OnePeriod until (nworms<5) or (t+tiny>=stoptime) or (overflow);
          CalculateMeanAndStdev;
          writeln(pherofile,round(t),',',dauerawakenprob:8:3,',',lambda:8:3,',',0.1*round(lambda*10):5:1,',',dauercount,',',mean:8:1,',',stdev:8:2);
*)
          for d:=0 to 4 do begin
          dauerawakenprob:=0.05-0.01*d;
          for ilambda:=0 to 50 do begin
            sumt:=0;
            writeln(d,ilambda:4, '/50');
            sumt:=0;
            for k:=1 to 30 do for c:=0 to 9 do begin
              feedfilename:='Ever_Sparser'+char(48+c)+'.inp';
              lambda:=-ilambda/12 + 0.5;
              adjustedfoodthreshold:=dauerfoodthresh;
              Init;
              if (woutput) then OpenOutputFile(1,1,1,1);
              repeat OnePeriod until (nworms<5) or (t+tiny>=stoptime) or (overflow);
              CalculateMeanAndStdev;
              inc(sumt,round(t));
              end;
            writeln(pherofile,sumt,',',lambda:8:3,',',round(t),',',dauerawakenprob:8:3,',',0.1*round(lambda*10):5:1,',',dauercount,',',mean:8:1,',',stdev:8:2);
            end;end;

          close(pherofile);
          end;
{$endif PHEROMONE}

{Main}
begin
{$ifdef PHEROMONE}
ReadInputFile;
OneMuRun(2,0.3);
OneMuRun(1,0.3);
OneMuRun(0,0.3);
OneMuRun(-1,0.3);
OneMuRun(-2,0.3);

OneMuRun(2,0.13);
OneMuRun(1,0.13);
OneMuRun(0,0.13);
OneMuRun(-1,0.13);
OneMuRun(-2,0.13);

OneMuRun(2,0.07);
OneMuRun(1,0.07);
OneMuRun(0,0.07);
OneMuRun(-1,0.07);
OneMuRun(-2,0.07);

(*
writeln(pherofile,'Varying Dauer-Awaken-Prob');
writeln(pherofile,crlf,'DauerAwakenProb,end-t, dauercount'{,' meanpop, stdev-pop'});
for ilambda:=1 to 200 do begin
  dauerawakenprob:=0.015+ilambda/1000;
  Init;
  {OpenOutputFile(1,1,1,1);}
  repeat OnePeriod until (nworms<5) or (t+tiny>=stoptime) or (overflow);
  {CalculateMeanAndStdev;}
  writeln(pherofile,dauerawakenprob:8:3,',',round(t),',',dauercount{,',',mean:8:1,',',stdev:8:2});
  end;
*)
end.
{$else PHEROMONE}
ReadInputFile;
TestFertility;end.
MainLoop;
write(allworms,',',matureworms,',',oldageworms);  Sleep(500);
end.
{$endif PHEROMONE}


(*  2017Sept08
procedure OnePeriod;
          type intarray=array of integer;
          var i,j,nworms_   :integer;
              {randorder     :^intarray;}
              interval      :double;

   {procedure Shuffle;
             var i,k,m :integer;
             begin
             for i:=1 to nworms_ do randorder^[i]:=i;
             for i:=1 to nworms_ do begin
                k:=randorder^[i];
                m:=succ(random(nworms_));
                randorder^[i]:=randorder^[m];
                randorder^[m]:=k;
                end;
             end;}

   procedure Shuffle;
             var i,m     :integer;
                 w1,w2   :worm;
             begin
             for i:=1 to nworms_ do begin
                w1:=w[i];
                m:=succ(random(nworms_));
                w[i]:=w[m];
                w[m]:=w1;
                end;
             end;

   procedure SumAllAppetites;
             var i :integer;
             begin
             if (foodadded<0) then food:=-foodadded;
             appsum:=0;
             for i:=1 to nworms do appsum:=appsum+w[i].appetite;
             end;


          begin {OnePeriod}
          if (feedinginterval=-1) then begin
            if (round(t)>=nextfeeding) then begin
              if (feedandculltogether) then FeedAndCull else Feed;
              ReadFromSchedFile;
              end;
            end
          else if (feedandculltogether) then begin
            if (feedinginterval=7223) or (feedinginterval=7322) then
              if ( (round(lastfeeding) mod 56 = 0) and (round(t-lastfeeding)>=24) )
                                          or
                 ( (round(lastfeeding) mod 56 > 0) and (round(t-lastfeeding)>=16) )
                 then FeedAndCull;
            if (round(t-lastfeeding)>=feedinginterval) then FeedAndCull; {if feedinginterval=7223 we'll never get here}
            end
          else begin
            if (cullinginterval=7223) or (cullinginterval=7322) then
              if ( (round(lastculling) mod 56 = 0) and (round(t-lastculling)>=24) )
                                          or
                 ( (round(lastculling) mod 56 > 0) and (round(t-lastculling)>=16) )
                 then Cull;
            if (round(t-lastculling)>=cullinginterval) then Cull; {if cullinginterval=7223 we'll never get here}
            if (feedinginterval=7223) or (feedinginterval=7322) then
              if ( (round(lastfeeding) mod 56 = 0) and (round(t-lastfeeding)>=24) )
                                          or
                 ( (round(lastfeeding) mod 56 > 0) and (round(t-lastfeeding)>=16) )
                 then Feed;
            if (round(t-lastfeeding)>=feedinginterval) then Feed; {if feedinginterval=7223 we'll never get here}
            end;
          DeListDeadWorms;
          if (nworms=0) then exit;
          nworms_:=nworms;
          {GetMem(randorder,sizeof(integer)*(nworms_+2));}
          interval:=one/nworms;
          {Shuffle;}
          SumAllAppetites;
          if (round(t)=100) and (sampletype<>'E') then samplecount:=0 {initializes record-keeping} else if (round(t)=200) then sampleon:=true; {initializes record-reporting}
          for i:=1 to nworms_ do begin
            t:=t+interval;
            {j:=randorder^[i];} j:=i;
            inc(w[j].age);
            if (foodadded<0) then food:=-foodadded;
            with w[j] do case s of  {each of these procedures has probabilities built in}
               egg    : Hatch;
               larva  : begin if Punt then continue; Eat; Alchemy; if not Mature then Starve; end;
               dauer  : if not Attrition then Awaken;
               adult  : begin Eat; Alchemy; Reproduce; if (not Bag) then SenescentMortality; end;
               wormbag: Burst;
               dead   : begin writeln('Dead worm found where there should not be one.'); readln; end;
               end;
            {$ifdef RDEBUG} if (t>100) and (t<110) then begin writeln(t:6:3,j:8); Checkup; end; {$endif RDEBUG}
            w[j].lastfood:=food;
            appsum:=appsum-w[j].appetite; if (appsum<one) then appsum:=one;
{$ifdef INDOUTPUT}
            with w[j] do if (samplecount<samplen) and (age<=maxlife) then RecordLHState;
{$endif INDOUTPUT}
            end;
          {FreeMem(randorder,sizeof(integer)*(nworms_+2));}
          DeListDeadWorms; {invokes Eliminate, where sample output is written to disk}
          RateHike(3,food*(1-foodattrition),0);
          food:=foodattrition*food;  {food loss, according to an input variable}
          {ComputeLarvaMass;
          if (mlarvae-oldmlarvae -(massr[4]+massr[13]-massr[9]+massr[6]-massr[7]-massr[8]-massr[10]-massr[11])>1E-5) then begin write(' End of Time Step'); readln; end;}
          ComputeStats;
          end;  {OnePeriod}
*)
