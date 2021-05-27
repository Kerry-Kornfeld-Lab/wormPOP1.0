unit DOSSTUFF;
{$J+}
{$define LAZARUS}
INTERFACE
uses SYSUTILS;

const breakkeys=[#3,#27,'Q','q']; whitespace=[#1..#32,'='];  crlf=#13#10;   numchars=['0'..'9','.','E','e','-'];
      Sunday=1; Monday=2; Tuesday=3; Wednesday=4; Thursday=5; Friday=6; Saturday=7;
      weekday:byte=0;
      thismonth:byte=0;
      day_abbrev:array[1..7] of string=('Sun','Mon','Tue','Wed','Thu','Fri','Sat');

type NoParams=procedure;

function Trim(instr :string):string;
function ToUpper(s :string):string;
procedure TrimToUpper(var s :string);
procedure Delay(msec :word);
procedure Readkey(var ki :char);
procedure Pause;
function TwoStr(x :word; zero :boolean):string;
function TimeStr:string;
procedure clrscr;
function itoa(x :longint):string;
function ftoaNoE(x :double; n,m :byte):string;
function ftoa(x :double; n :byte):string;
function PadItoa(x :longint; len :byte):string;
function RightPadItoa(x :longint; len :byte):string;
procedure ClrEOL;
procedure LineReturn;
procedure NextFileName(var s :string);
procedure ActiveDelay(m :word; proc :NoParams); {m is time in milliseconds; proc is what to do while waiting}
function floor(x :double):integer;

IMPLEMENTATION

{$ifdef LAZARUS}
procedure ActiveDelay(m :word; proc :NoParams); {m is time in milliseconds; proc is what to do while waiting}
          const ms=1.0/24/3600/1000;
          var til :TDateTime;
          begin
          Sleep(m);
          end;
{$else LAZARUS}
procedure ActiveDelay(m :word; proc :NoParams); {m is time in milliseconds; proc is what to do while waiting}
          const ms=1.0/24/3600/1000;
          var til :TDateTime;
          begin
          til:=GetTime+m*ms;
          repeat proc until (GetTime=til);
          end;
{$endif LAZARUS}

procedure Pause;
          var ki :char;
          begin
          read(ki); if (ki in breakkeys) then halt;
          end;

procedure ReadKey(var ki :char);
          begin
          read(ki);
          end;

function Trim(instr :string):string;
         var i :byte;
         begin
         result:=instr;
         i:=0;
         while (result[succ(i)]) in whitespace do inc(i);
         if (i>0) then Delete(result,1,i);
         i:=succ(length(result));
         while (result[pred(i)]) in whitespace do dec(i);
         if (i<=length(result)) then Delete(result,i,succ(length(result)-i));
         end;

function ToUpper(s :string):string;
         var i :byte;
         begin
         result:='';
         for i:=1 to length(s) do result:=result+upcase(s[i]);
         end;

function floor(x :double):integer;
         begin
         result:=round(x-0.5);
         end;

procedure TrimToUpper(var s :string);
          var i :byte;
          begin
          i:=0;
          while (s[succ(i)]) in whitespace do inc(i);
          if (i>0) then Delete(s,1,i);
          i:=succ(length(s));
          while (s[pred(i)]) in whitespace do dec(i);
          if (i<=length(s)) then Delete(s,i,succ(length(s)-i));
          for i:=1 to length(s) do s[i]:=upcase(s[i]);
          end;

function ftoa(x :double; n :byte):string; {n is specified total string length}
         var p      :byte;
         begin
         if (n<4) then n:=4;
         Str(x:n+6,result);
         p:=2+Pos('E',result);
         while (result[p]='0') and (p<length(result)) do Delete(result,p,1);
         if result[p-1]='+' then Delete(result,p-1,1);
         dec(p,2);
         if (length(result)>n) then repeat
           Delete(result,p-1,1);
           p:=Pos('E',result);
         until (length(result)=n) or (result[p-1]='.');
         end;

function ftoaNoE(x :double; n,m :byte):string; {n is specified total string length, m is #decimals}
         begin
         if (x>1E10) then begin result:='Infinity'; if (n<8) then result:=copy(result,1,n); exit; end; 
         Str(x:n+6:m,result);
         while (result[1]=' ') and (length(result)>n) do delete(result,1,1);
         if (n=0) then exit; {n=0 means indeterminate length}
         while (result[length(result)]>'.') and (length(result)>n) do SetLength(result,pred(length(result)));
         if (length(result)>n) and (result[length(result)]='.') then SetLength(result,pred(length(result)));
         end;

function itoa(x :longint):string;
         var i    :integer;
         begin
         Str(x,result);
         i:=length(result)-2;
         while (i>1) do begin
           Insert(',',result,i);
           dec(i,3);
           end;
         itoa:=result;
         end;

function PadItoa(x :longint; len :byte):string;
         begin
         result:=itoa(x);
         while (length(result)<len) do result:=' '+result;
         end;

function RightPadItoa(x :longint; len :byte):string;
         begin
         result:=itoa(x);
         while (length(result)<len) do result:=result+' ';
         end;

procedure Delay(msec :word);
          const daysperms=1/(24*3600*1000);
          var timenow, targettime :TDateTime;
          begin
          targettime:=Time+msec*daysperms;
          repeat
            timenow:=time;
          until (timenow>targettime);
          end;

function TwoStr(x :word; zero :boolean):string;
         var ws :string;
         begin
         Str(x:2,ws);
         if (ws[1]=' ') then
           if (zero) then ws[1]:='0' else ws:=ws[2];
         TwoStr:=ws;
         end;

function TimeStr:string;
         begin
         TimeStr:=DateTimeToStr(date+time);
         end;
(*
function TimeStr:string;
         var days :TDateTime;
         begin
         days:=Time;
         TimeStr:=FormatDateTime('c',days);
         end;
*)
procedure clrscr;
          var i :integer;
          begin
          for i:=1 to 25 do writeln;
          end;

procedure ClrEOL;
        begin
        write(StringOfChar(' ',79));
        write(StringOfChar(#8,79));
        end;

procedure LineReturn;
        begin
        write(#13); ClrEOL;
        end;

procedure IncrementChar(var s: string; p :byte);
          var ch :char;
              next :boolean;
          begin
          ch:=s[p];
          next:=true;
          if (ch='z') then ch:='a'
          else if (ch='Z') then ch:='A'
          else if (ch='9') then ch:='0'
          else begin next:=false; inc(ch); end;
          {s:=StuffString(s,p,1,ch);}
          s[p]:=ch;
          if (next) then IncrementChar(s,pred(p));
          end;

procedure NextFileName(var s :string);
          begin
          while FileExists(s) do IncrementChar(s,length(s));
          end;

var ta:integer;

begin
weekday:=DayOfWeek(Date);
val(FormatDateTime('MM',date-8),thismonth,ta);
end.
