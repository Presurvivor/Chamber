clear variables
close all
clc
tic;
%% FILE NAME DECLARED IN CODE ---------------------------------------------

fileToRead = 'Book1.xlsx';
% fileToRead = 'System-Check.IN617.v1.xlsx';

% FILE NAME READ FROM USER INPUT------------------------------------------
%fileToRead = input('Enter the file name (including the extension):\n', 's');
% -------------------------------------------------------------------------

InputFile = fileToRead(1:(length(fileToRead)-5));

newFolder = strcat(InputFile, ' RESULTS');

mkdir(newFolder);

imp_data = xlsread(fileToRead);

    %Checking for non-inclution of headers
    if isnan(imp_data(1,1)) == 1 && isnan(imp_data(2,1))== 1
        imp_data(1:2,:) = [];
    end
% -------------------------------------------------------------------------

case_number_A = imp_data(:,1);
method_A = imp_data(:,2);
write_A = imp_data(:,3);
mat_ID_A = imp_data(:,4);

%Mechanical inputs
MinMech_A = imp_data(:,6);
MaxMech_A = imp_data(:,7);
RateMech_A = imp_data(:,8);
Hold_Tens_A = imp_data(:,9);
Hold_Comp_A = imp_data(:,10);
Soak_time_A = imp_data(:,11);

%Temperature inputs
MaxT_A = imp_data(:,12);
MinT_A = imp_data(:,13);
AmbientT_A = imp_data(:,14);
Phase_Lag_A = imp_data(:,15);
Phasing_A = imp_data(:,16);

% Acoustical Inputs
SPL_a_A = imp_data(:,17);
SPL_Dwell_A = imp_data(:,18);
SPL_f_A = imp_data(:,19);
Phase_Lag_Acoustic_A = imp_data(:,20);

t_inc_A = imp_data(:,21);
num_cycles_A = imp_data(:,22);

pi= 3.14159;

for i=1:length(case_number_A);
%% Mechanical Waveform ========================================================================================

case_number = case_number_A(i);
method = method_A(i);
write = write_A(i);
mat_ID = mat_ID_A(i);
MinMech = MinMech_A(i);
MaxMech = MaxMech_A(i);
RateMech = RateMech_A(i);
Hold_Tens = Hold_Tens_A(i);
Hold_Comp = Hold_Comp_A(i);
Soak_time = Soak_time_A(i);
MaxT = MaxT_A(i);
MinT = MinT_A(i);
AmbientT = AmbientT_A(i);
Phase_Lag = Phase_Lag_A(i);
Phasing = Phasing_A(i);
SPL_a = SPL_a_A(i);
SPL_Dwell = SPL_Dwell_A(i);
SPL_f = SPL_f_A(i);
Phase_Lag_Acoustic = Phase_Lag_Acoustic_A(i);
t_inc = t_inc_A(i);
num_cycles = num_cycles_A(i);


Ratepercentcps = RateMech*100;
Ratepercentpmin = Ratepercentcps*60;

MeanMech = 0.5*(MaxMech+MinMech);
AmpMech = 0.5*(MaxMech-MinMech);
RangeMech = MaxMech-MinMech;

R = MinMech/(MaxMech+0.00000001);

AMech = (MaxMech-MinMech)/(0.000000001+MaxMech+MinMech);
MMech = (MaxMech+MinMech)/(MaxMech-MinMech);
t_cyc = ((2*RangeMech)/(RateMech))+Hold_Tens+Hold_Comp;

t = 0:t_inc:(t_cyc*num_cycles);

            %REDUCING t IN TENSILE DWELL PERIOD====================================

            if (Hold_Tens ~= 0)

                t_esc = 1.05;
                T = t;
                for n = 1:num_cycles;

                    Tti = Soak_time + (t_cyc)*(n-1) + (MaxMech/RateMech);
                    Ttf = Tti + Hold_Tens;

                    TDwell_length = Tti + t_inc;

                    d_t_inc = t_inc;        %WHEN ON, ..... . . .  ||  ..... . . .

                    while (TDwell_length(end) < Ttf)
                        d_t_inc = t_esc * d_t_inc;
                        TDwell_length(end + 1) = TDwell_length(end) + d_t_inc;                
                    end


                    TDwell_length(end) = Ttf;
                    Tti_index = find(T == Tti);
                    Ttf_index = find(T == Ttf);

                    T = [T(1:Tti_index) TDwell_length T((Ttf_index + 1):end)];
                end
                t = T;
            end

            %======================================================================

            % REDUCING t IN COMPRESSIVE DWELL PERIOD ==============================

            if (Hold_Comp ~= 0)
                t_esc = 1.05;
                T = t;
                for n = 1:num_cycles;

                    Cti = Soak_time + (t_cyc)*(n-1) + (MaxMech/RateMech) + Hold_Tens + (RangeMech/RateMech);
                    Ctf = Cti + Hold_Comp;

                    CDwell_length = Cti + t_inc;

                    d_t_inc = t_inc;        %WHEN ON, ..... . . .  ||  ..... . . .

                    while (CDwell_length(end) < Ctf)
                        d_t_inc = t_esc * d_t_inc;
                        CDwell_length(end + 1) = CDwell_length(end) + d_t_inc;                
                    end

                    CDwell_length(end) = Ctf;
                    Cti_index = find(abs(T-Cti)<t_inc);
                    if (length(Cti_index) == 2)
                        Cti_index = Cti_index(2);
                    end
                    
                    Ctf_index = find(abs(T-Ctf)<t_inc);
                    if (length(Ctf_index) == 2)
                        Ctf_index = Ctf_index(2);
                    end

                    T = [T(1:Cti_index) CDwell_length T((Ctf_index + 1):end)];
                end
                t = T;

            end
            %======================================================================    

%Dwell Fatigue
Del_Comp = -RateMech*((Hold_Comp)/2);
Del_Tens = RateMech*((Hold_Tens)/2);
AmpMech2 = (MaxMech+Del_Tens-MinMech-Del_Comp)/2;
MeanMech2 = (MaxMech+Del_Tens+MinMech+Del_Comp)/2;

%Functional Description
WaitMech = sign(MaxMech+MinMech+0.00000001)*min(((MeanMech2)/RateMech),((abs(MeanMech2))/RateMech));
WaitMech2 = (abs(MaxMech-MinMech))/(2*RateMech);


Mechanical1 = max(min(MaxMech,MeanMech2+((2*AmpMech2*asin((sin(pi*t*(RateMech*sign(MeanMech+0.000000001)/(2*AmpMech2))))))/pi)),MinMech);

Mechanical2 = (RateMech*t*sign(MeanMech+0.00000000001));

if t > WaitMech2
    Mechanical3 = max(min(MaxMech,MeanMech2+((2*AmpMech2*asin((sin(pi*(t-WaitMech)*(RateMech*sign(MeanMech+0.000000001)/(2*AmpMech2))))))/pi)),MinMech);
else
    Mechanical3 = (RateMech*t*sign(MeanMech+0.00000000001));
end

if t < Soak_time                                                                        %  WHY IS THIS NEEDED?
            Mechanical4 = 0;                                                            %
elseif t < WaitMech2                                                                    %
            Mechanical4 = (RateMech*(t-Soak_time)*sign(MeanMech+0.00000000001));        %
else
            Mechanical4 = max(min(MaxMech,MeanMech2+((2*AmpMech2*asin((sin(pi*(t-Soak_time-WaitMech)*(RateMech*sign(MeanMech+0.000000001)/(2*AmpMech2))))))/pi)),MinMech);
end

%CHANGES MADE ON 3/10/2015 ------------------------------------------------

% PHASING CHANGE 
RangeT = MaxT - MinT;

if (Phasing == 180 || Phasing == -180) && (RangeT ~= 0) && (MinMech == -MaxMech)
    Mechanical4 = -Mechanical4;
end


%SUPRESSING A PART OF WAVEFORM
 
check1_t = Soak_time + (2*t_inc);
check2_t = Soak_time + (3*t_inc);

check1_t_index = ind2sub(size(t),find(t == check1_t));
check2_t_index = ind2sub(size(t),find(t == check2_t));

if (Hold_Tens == 0) && (Hold_Comp ~= 0) && (Mechanical4(check1_t_index) < Mechanical4(check2_t_index))
    
    aCheck = Soak_time + (MaxMech/RateMech) + ((RangeMech/RateMech)/2);
    aCheck_index = find(abs(t-aCheck)<t_inc);
    aCheck_index = aCheck_index(2);
    Mechanical4(2:aCheck_index) = [];
    subtractVal = t(aCheck_index);
    t(1:aCheck_index-1) = [];
    t = t - subtractVal;
    
end

%--------------------------------------------------------------------------

%% Temperature Waveform =================================================================================================================

Phase_time = (Phase_Lag)/(2*pi)*(t_cyc);
% RangeT = MaxT - MinT;
MeanT = (.5)*(MaxT+MinT);

if t_cyc == 0
RateT = RateTemp;
else
RateT =(2)*(RangeT/t_cyc);
end

Del_Comp_T = (-RateT)*((Hold_Comp)/2);
Del_Tens_T = (RateT)*((Hold_Tens)/2);
Amp2T = (((MaxT)+(Del_Tens_T) - (MinT)-(Del_Comp_T))/2);
Mean2T = (((MaxT)+(Del_Tens_T)+(MinT)+(Del_Comp_T))/2);

if t_cyc == 0
PeriodT = ((2*RangeT)/(RateTemp))+Hold_Tens+Hold_Comp;
else
PeriodT = t_cyc;        %  REDUNDANT BC t = 0:t_inc:(t_cyc*num_cycles);
end                     
                        %  ...... THEREFORE......

                        % t = 0:t_inc:(PeriodT*num_cycles);


%fix Mechanical3(t-Phase_time)
if t > WaitMech2
            Mechanical3_1 = max(min(MaxMech,MeanMech2+((2*AmpMech2*asin((sin(pi*((t-Phase_time)-WaitMech)*(RateMech*sign(MeanMech+0.000000001)/(2*AmpMech2))))))/pi)),MinMech);
else
            Mechanical3_1 = (RateMech*(t-Phase_time)*sign(MeanMech+0.00000000001));
end
%waveform continued
Temp1 = (((Mechanical3_1)-MinMech)/(MaxMech-MinMech))*(MaxT-MinT)*((2*Phasing-pi)/(-pi))+((MinT)*((pi-Phasing)/pi))+((MaxT)*(Phasing/pi)); 
%fix Temp1(Phase_time)
if t > WaitMech2
            Mechanical3_2 = max(min(MaxMech,MeanMech2+((2*AmpMech2*asin((sin(pi*((Phase_time-Phase_time)-WaitMech)*(RateMech*sign(MeanMech+0.000000001)/(2*AmpMech2))))))/pi)),MinMech);
else
            Mechanical3_2 = (RateMech*(Phase_time-Phase_time)*sign(MeanMech+0.00000000001));
end
        
%waveform continued
if	t < Phase_time
		Temp2 =  (((Mechanical3_2)-MinMech)/(MaxMech-MinMech))*(MaxT-MinT)*((2*Phasing-pi)/(-pi))+((MinT)*((pi-Phasing)/pi))+((MaxT)*(Phasing/pi)); 
else	
		Temp2 = (((Mechanical3_1)-MinMech)/(MaxMech-MinMech))*(MaxT-MinT)*((2*Phasing-pi)/(-pi))+((MinT)*((pi-Phasing)/pi))+((MaxT)*(Phasing/pi)); 
end
	
%fix Mechanical3(t-Phase_time-PeriodT/4)
if t > WaitMech2
            Mechanical3_3 = max(min(MaxMech,MeanMech2+((2*AmpMech2*asin((sin(pi*((((t)-(Phase_time)-(PeriodT/4))-Phase_time)-WaitMech)*(RateMech*sign(MeanMech+0.000000001)/(2*AmpMech2))))))/pi)),MinMech);
else
            Mechanical3_3 = (RateMech*(((t)-(Phase_time)-(PeriodT/4))-Phase_time)*sign(MeanMech+0.00000000001));
end

%waveform continues    
Temp3 = ((((Mechanical3_3)-(MinMech))/(MaxMech - MinMech))*((MaxT-MinT)+(MinT)));

%fix mechanical3(Phase_time + (PeriodT/4))

if t > WaitMech2
            Mechanical3_4 = max(min(MaxMech,MeanMech2+((2*AmpMech2*asin((sin(pi*(((Phase_time + (PeriodT/4))-Phase_time)-WaitMech)*(RateMech*sign(MeanMech+0.000000001)/(2*AmpMech2))))))/pi)),MinMech);
else
            Mechanical3_4 = (RateMech*((Phase_time + (PeriodT/4))-Phase_time)*sign(MeanMech+0.00000000001));
end

%waveform continues
if	t < (Phase_time + (PeriodT/4))
		Temp4 =(((((Mechanical3_4)-(Phase_time)-(PeriodT/4))-(MinMech))/(MaxMech - MinMech))*((MaxT-MinT)+(MinT)));
	
else
		Temp4 = ((((Mechanical3_4)-(MinMech))/(MaxMech - MinMech))*((MaxT-MinT)+(MinT)));
end	
 
%Temp4(0)= Temp4_0=Temp4   
   
if Phasing == 90;          % (3/10/2015) SHOULD (3.14159/2) BE CHANGED TO 90?
		Temp5 = Temp4;
elseif	t < Phase_time
		Temp5 =  (((Mechanical3_2)-MinMech)/(MaxMech-MinMech))*(MaxT-MinT)*((2*Phasing-pi)/(-pi))+((MinT)*((pi-Phasing)/pi))+((MaxT)*(Phasing/pi)); 
else	
		Temp5 = (((Mechanical3_1)-MinMech)/(MaxMech-MinMech))*(MaxT-MinT)*((2*Phasing-pi)/(-pi))+((MinT)*((pi-Phasing)/pi))+((MaxT)*(Phasing/pi)); 
end
   		
Temp6 = max(min(MaxT,Mean2T+((2*Amp2T*asin((sin(pi*t*(RateT*sign(MeanT+0.000000001)/(2*Amp2T))))))/pi)),MinT);

if	MaxMech == 0 ^ MinMech == 0
		Temp7 = max(min(MaxT,Mean2T+((2*Amp2T*asin((sin(pi*t*(RateT*sign(MeanT+0.000000001)/(2*Amp2T))))))/pi)),MinT);
elseif Phasing == (3.14159/2);
		Temp7 = Temp4;
elseif	t < Phase_time
		Temp7 =  (((Mechanical3_2)-MinMech)/(MaxMech-MinMech))*(MaxT-MinT)*((2*Phasing-pi)/(-pi))+((MinT)*((pi-Phasing)/pi))+((MaxT)*(Phasing/pi)); 
else	
		Temp7 = (((Mechanical3_1)-MinMech)/(MaxMech-MinMech))*(MaxT-MinT)*((2*Phasing-pi)/(-pi))+((MinT)*((pi-Phasing)/pi))+((MaxT)*(Phasing/pi)); 
end
	
	
if Temp7 > MaxT;
		Temp8 = MaxT;
elseif	MaxMech == 0 ^ MinMech == 0
		Temp8 = max(min(MaxT,Mean2T+((2*Amp2T*asin((sin(pi*t*(RateT*sign(MeanT+0.000000001)/(2*Amp2T))))))/pi)),MinT);
elseif Phasing == (3.14159/2);
		Temp8 = Temp4;
elseif	t < Phase_time;
		Temp8 =  (((Mechanical3_2)-MinMech)/(MaxMech-MinMech))*(MaxT-MinT)*((2*Phasing-pi)/(-pi))+((MinT)*((pi-Phasing)/pi))+((MaxT)*(Phasing/pi)); 
else	
		Temp8 = (((Mechanical3_1)-MinMech)/(MaxMech-MinMech))*(MaxT-MinT)*((2*Phasing-pi)/(-pi))+((MinT)*((pi-Phasing)/pi))+((MaxT)*(Phasing/pi)); 
end
    

if Temp8 < MinT
		Temp9 = MinT;
elseif Temp7 > MaxT;
		Temp9 = MaxT;
elseif	MaxMech == 0 ^ MinMech == 0
		Temp9 = max(min(MaxT,Mean2T+((2*Amp2T*asin((sin(pi*t*(RateT*sign(MeanT+0.000000001)/(2*Amp2T))))))/pi)),MinT);
elseif Phasing == (3.14159/2);
		Temp9 = Temp4;
elseif	t < Phase_time
		Temp9 =  (((Mechanical3_2)-MinMech)/(MaxMech-MinMech))*(MaxT-MinT)*((2*Phasing-pi)/(-pi))+((MinT)*((pi-Phasing)/pi))+((MaxT)*(Phasing/pi)); 
else	
		Temp9 = (((Mechanical3_1)-MinMech)/(MaxMech-MinMech))*(MaxT-MinT)*((2*Phasing-pi)/(-pi))+((MinT)*((pi-Phasing)/pi))+((MaxT)*(Phasing/pi)); 
end
      
 %Mechanical3_1(0)
 if t > WaitMech2
            Mechanical3_1 = max(min(MaxMech,MeanMech2+((2*AmpMech2*asin((sin(pi*((0-Phase_time)-WaitMech)*(RateMech*sign(MeanMech+0.000000001)/(2*AmpMech2))))))/pi)),MinMech);
 else
            Mechanical3_1 = (RateMech*(0-Phase_time)*sign(MeanMech+0.00000000001));
 end
 %Mechanical3_1(t-soaktime)
 if t > WaitMech2
            Mechanical3_1_Soak_time = max(min(MaxMech,MeanMech2+((2*AmpMech2*asin((sin(pi*(((t-Soak_time)-Phase_time)-WaitMech)*(RateMech*sign(MeanMech+0.000000001)/(2*AmpMech2))))))/pi)),MinMech);
 else
            Mechanical3_1_Soak_time = (RateMech*((t-Soak_time)-Phase_time)*sign(MeanMech+0.00000000001));
 end
 %Mechanical3_2(0)= Mechanical3_t(t-soaktime)= Mechanical3_2
 
 %Temp9(0)= Temp9_0
 
if Temp8 < MinT;
		Temp9_0 = MinT;
elseif Temp7 > MaxT
		Temp9_0 = MaxT;
elseif	MaxMech == 0 ^ MinMech == 0
		Temp9_0 = max(min(MaxT,Mean2T+((2*Amp2T*asin((sin(pi*t*(RateT*sign(MeanT+0.000000001)/(2*Amp2T))))))/pi)),MinT);
elseif Phasing == (3.14159/2);
		Temp9_0 = Temp4;
elseif	t < Phase_time
		Temp9_0 =  (((Mechanical3_2)-MinMech)/(MaxMech-MinMech))*(MaxT-MinT)*((2*Phasing-pi)/(-pi))+((MinT)*((pi-Phasing)/pi))+((MaxT)*(Phasing/pi)); 
else	
		Temp9_0 = (((Mechanical3_1_Soak_time)-MinMech)/(MaxMech-MinMech))*(MaxT-MinT)*((2*Phasing-pi)/(-pi))+((MinT)*((pi-Phasing)/pi))+((MaxT)*(Phasing/pi)); 
end
        
    %Temp9(t-Soak_time)
if Temp8 < MinT
		Temp9_Soak_time = MinT;
elseif Temp7 > MaxT
		Temp9_Soak_time = MaxT;
elseif	MaxMech == 0 ^ MinMech == 0
		Temp9_Soak_time = max(min(MaxT,Mean2T+((2*Amp2T*asin((sin(pi*t*(RateT*sign(MeanT+0.000000001)/(2*Amp2T))))))/pi)),MinT);
elseif Phasing == (3.14159/2);
		Temp9_Soak_time = Temp4;
elseif	t < Phase_time
		Temp9_Soak_time =  (((Mechanical3_2)-MinMech)/(MaxMech-MinMech))*(MaxT-MinT)*((2*Phasing-pi)/(-pi))+((MinT)*((pi-Phasing)/pi))+((MaxT)*(Phasing/pi)); 
else	
		Temp9_Soak_time = (((Mechanical3_1_Soak_time)-MinMech)/(MaxMech-MinMech))*(MaxT-MinT)*((2*Phasing-pi)/(-pi))+((MinT)*((pi-Phasing)/pi))+((MaxT)*(Phasing/pi)); 
end

%Changed Temp9_Soak_time to Temp9_Soak_time
 %Waveform continued
if t < Soak_time
        Temp10 = Temp9_0;
else
        Temp10 = Temp9_Soak_time;
end


vconvert = Temp10*((10+0)/((1368)-(-270)))+1.65;

%% Acoustical Waveform ========================================================================================

if t_cyc == 0
    PeriodA = 1;
else
    PeriodA = t_cyc;
end


Phase_Time_Acoustic = (((Phase_Lag_Acoustic)/(2*pi))*(t_cyc));

%Mechanical4(t+t_inc-Phase_Time_Acoustic)
if t < Soak_time
            Mechanical41 = 0;
elseif t > WaitMech2
            Mechanical41= max(min(MaxMech,MeanMech2+((2*AmpMech2*asin((sin(pi*(((t+t_inc-Phase_Time_Acoustic))-Soak_time-WaitMech)*(RateMech*sign(MeanMech+0.000000001)/(2*AmpMech2))))))/pi)),MinMech);
else
            Mechanical41 = (RateMech*(((t+t_inc-Phase_Time_Acoustic))-Soak_time)*sign(MeanMech+0.00000000001));
end
%Mechanical4(t-Phase_Time_Acoustic)
if t < Soak_time
            Mechanical42= 0;
elseif t > WaitMech2
            Mechanical42= max(min(MaxMech,MeanMech2+((2*AmpMech2*asin((sin(pi*((t-Phase_Time_Acoustic)-Soak_time-WaitMech)*(RateMech*sign(MeanMech+0.000000001)/(2*AmpMech2))))))/pi)),MinMech);
else
            Mechanical42 = (RateMech*((t-Phase_Time_Acoustic)-Soak_time)*sign(MeanMech+0.00000000001));
end

 %Acoustical Waveform Continued
 
 SPL_t = SPL_a*(1-SPL_Dwell*(abs(sign((Mechanical41)-(Mechanical42)))));

 %Mechanical41(t-Soak_time) = Mechanical43
 if t < Soak_time
            Mechanical43 = 0;
elseif t > WaitMech2
            Mechanical43= max(min(MaxMech,MeanMech2+((2*AmpMech2*asin((sin(pi*((((t-Soak_time)+t_inc-Phase_Time_Acoustic))-Soak_time-WaitMech)*(RateMech*sign(MeanMech+0.000000001)/(2*AmpMech2))))))/pi)),MinMech);
else
            Mechanical43 = (RateMech*((((t-Soak_time)+t_inc-Phase_Time_Acoustic))-Soak_time)*sign(MeanMech+0.00000000001));
end
 %Mechanical42(t-Soak_time) = Mechanical44
 if t < Soak_time
            Mechanical44 = 0;
elseif t > WaitMech2
            Mechanical44= max(min(MaxMech,MeanMech2+((2*AmpMech2*asin((sin(pi*(((t-Soak_time)-Phase_Time_Acoustic)-Soak_time-WaitMech)*(RateMech*sign(MeanMech+0.000000001)/(2*AmpMech2))))))/pi)),MinMech);
else
            Mechanical44 = (RateMech*(((t-Soak_time)-Phase_Time_Acoustic)-Soak_time)*sign(MeanMech+0.00000000001));
 end

 %Acoustical Waveform Continued
 if t < Soak_time
     SPL_t2 = 0;
 else
     SPL_t2 = SPL_a*(1-SPL_Dwell*(abs(sign((Mechanical43)-(Mechanical44)))));
 end

 %% OUTPUTTING Waveform ========================================================================================

plot_name = sprintf('%d-WaveformPLOT',case_number);
f0 = fullfile(newFolder, plot_name);
h=figure('visible','off');
plot(t,Mechanical4);
print(f0,'-djpeg');

wfoutput = zeros(length(t),5);
wfoutput(:,1) = t;
wfoutput(:,2) = Mechanical4;
wfoutput(:,3) = Temp10;
wfoutput(:,4) = SPL_a;
wfoutput(:,5) = SPL_f;
wfoutput(:,6) = RateT;
 if write == 1 
     
     file_waveform = sprintf('%d-Waveform.csv', case_number);

     f1 = fullfile(newFolder, file_waveform);
     
     header_waveform = {'Time t (sec)' ,'Strain (mm/mm)','Temperature T (degC)', 'Vibration Amplitude (mm/mm)', 'Frequency f (Hz)'};
     
     csvwriteh(f1, wfoutput(:,1:5), header_waveform)
 
 end
 
end
toc;
