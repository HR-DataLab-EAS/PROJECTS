%%
% ##### https://nl.mathworks.com/help/bioinfo/ug/preprocessing-raw-mass-spectrometry-data.html
% ##### https://nl.mathworks.com/help/bioinfo/ug/visualizing-and-preprocessing-hyphenated-mass-spectrometry-data-sets-for-metabolite-and-protein-peptide-profiling.html
clc
clear all
% dir_data = 'D:\RAW\EXAMPLE\20241022';

dir_data = 'D:\RAW\EXAMPLE\20241022_TEST';
cd(dir_data)

% =====>  INLEZEN .raw converted to .mzXML file data WITH MATLAB

% % is also possible using Python
% % ===> https://github.com/elnurgar/mzxml-precursor-corrector

% % MATLAB REQUIRES mzXML FORMAT
% % CONVERSION from RAW to  mzXML 
% % requires proteowizard to be installed (windows 11) from https://proteowizard.sourceforge.io/download.html
% % with MSConvertGUI.exe

%%%%% ===> see ALSO https://github.com/awbirdsall/analyzems


% % set current directory to the map were the mzXML files are located
%%%%%%%%% cd('D:\RAW\PROTEO_WIZZARD\rawdata')
cd(dir_data)

% % Matlab: use functions from the Bioinformatics Toolbox
%%% ====> https://nl.mathworks.com/help/bioinfo/ref/mzxml2peaks.html
%%% ====> https://nl.mathworks.com/help/bioinfo/ug/preprocessing-raw-mass-spectrometry-data.html
% % help mzxml2peaks
tic
% start timer

% ##### FILENAME #### % 
mzfile = 'Lemonfiness_s3_onbespoten.mzxml';

% ##### META INFO #### % 
mzinfo = mzxmlinfo(mzfile,'NUMOFLEVELS',true);

% Assuming mzinfo is the struct obtained from mzxmlinfo(mzfile)

% Extract field names and values
fieldNames = fieldnames(mzinfo);
numFields = length(fieldNames);
fieldValues = cell(numFields, 1);

for i = 1:numFields
    fieldName = fieldNames{i};
    fieldValues{i} = mzinfo.(fieldName); % Access field value using dot notation
end

% Create the table
Table_mzInfo = table(fieldNames, fieldValues, 'VariableNames', {'FieldName', 'FieldValue'});

% Display the table
disp(Table_mzInfo);

% ##### mzxml comtains 3 fields  scan | mzXML | index #### % 
mzxml_struct = mzxmlread(mzfile);

% ##### The ROOT element of our schema is called mzXML #### % 
mzXML = mzxml_struct.mzXML;

% ### msInstrument #### % 
% ===>  resolution, manufacturer, model, ionization type, 
% ===>   mass analyzer type, detector type) and acquisition software used to generate the data

field = fieldnames(mzXML.msRun.msInstrument);
numFields = length(field);

% Preallocate cell arrays for efficiency
fieldNames = cell(numFields, 1);
fieldValues = cell(numFields, 1);

        for i = 1:numFields
            fieldName = field{i};
            fieldValue = getfield(mzXML.msRun.msInstrument, fieldName);
        
            % Handle potential structures within fieldValue
            if isstruct(fieldValue)
                % Extract the 'value' field if it exists
                if isfield(fieldValue, 'value')
                    fieldValue = fieldValue.value;
                elseif strcmp(fieldName, 'software') % Handle software field specifically
                    fieldValue = [fieldValue.name, ' ', fieldValue.version]; 
                else
                    % If no 'value' field, handle differently (e.g., convert to string)
                    fieldValue = struct2table(fieldValue); 
                end
            end
        
            fieldNames{i} = fieldName;
            fieldValues{i} = fieldValue;
        end

% Create the table
Table_instrument = table(fieldNames, fieldValues, 'VariableNames', {'FieldName', 'FieldValue'});

% Display the table
disp(Table_instrument);


% % .scan contains 27 parameters
scan_data = mzxml_struct.scan;
params = fieldnames(mzxml_struct.scan);

% Create a table with item numbers and label names
Table_scan_params = table((1:numel(params))', params, 'VariableNames', {'ItemNumber', 'LabelName'});

% Display the table
disp(Table_scan_params)

% Get a list of all variables in the workspace
allVars = who;

% Find the indices of variables that start with "Table_"
tableVarIndices = find(startsWith(allVars, 'Table_'));

% Create a cell array with the names of the variables to keep
variablesToKeep = allVars(tableVarIndices);

% Clear all variables except those starting with "Table_"
clearvars('-except', variablesToKeep{:});

% ##### FILENAME #### % 
mzfile = 'Lemonfiness_s3_onbespoten.mzxml';

% stop timer + show the elapsed time
toc

%%
% ################################################################################## %
% ################################################################################## %
% ################################################################################## %
% ################################################################################## %
% ################################################################################## %
% ################################################################################## %
% ################################################################################## %
% ################################################################################## %
% ################################################################################## %
% ################################################################################## %
% ################################################################################## %

%%
% The |MZXMLREAD| function reads the XML document into a MATLAB structure.
% The fields |scan| and |index| are placed at the first level of the output
% structure for improved access to the spectral data. The remainder of the
% mzXML document tree is parsed according to the schema specifications.
% This LC/MS data set contains 528 +  4224 scans with two MS levels. 
% For this example you will use only the first level scans. 
% Second level spectra are usually used for peptide/protein identification, 
% and come at a later stage in some types of workflow analyses. 
% |MZXMLREAD| can filter the desired scans without loading all the dataset into memory:

clc
tic

% ##### FILENAME #### % 
cd('D:\RAW\EXAMPLE\20241022_TEST')
mzfile = 'Lemonfiness_s3_onbespoten.mzxml';

% ###################### MS1 scan ###################### %
% ###################### MS1 scan ###################### %
% ###################### MS1 scan ###################### %
% ###################### MS1 scan ###################### %
mzXML_struct_ms1 = mzxmlread(mzfile,'LEVEL',1);
ms1_scan_data = mzXML_struct_ms1.scan;
% ====> clear mzXML_struct_ms1

% Define the fields to keep
fieldsToKeep = {'num', 'msLevel', 'peaksCount', 'retentionTime', ...
               'startMz', 'endMz', 'lowMz', 'highMz', ...
               'basePeakMz', 'basePeakIntensity', 'totIonCurrent', 'peaks', 'msInstrumentID'};

% Remove all other fields from the struct array
ms1_scan_data = rmfield(ms1_scan_data, setdiff(fieldnames(ms1_scan_data), fieldsToKeep));
ms1_scan_xlsx = ms1_scan_data;
scan_num1 = cell2mat({ms1_scan_data.num})';

% ###################### MS1 peak data + retention time in seconds ###################### %
% ###################### MS1 peak data + retention time in seconds ###################### %
% ###################### MS1 peak data + retention time in seconds ###################### %
%%%%% TABELIZE PEAKS DATA  into [INTENSITY] absorbance [as a function of] m/z [MASS]
[ms1_peaks,ms1_ret_time] = mzxml2peaks(mzXML_struct_ms1, 'Levels', 1);

for i = 1:length(ms1_scan_data)
    %% DATA needed to create readable  .xlsx  files
    ms1_scan_xlsx(i).retentionTime = ms1_ret_time(i);

    %% DATA needed to create readable  .mat  files
    ms1_scan_data(i).retentionTime = ms1_ret_time(i); 
    ms1_scan_data(i).peaks = ms1_peaks(i); 
end

% ###################### MS2 scan ###################### %
% ###################### MS2 scan ###################### %
% ###################### MS2 scan ###################### %
% ###################### MS2 scan ###################### %
mzXML_struct_ms2 = mzxmlread(mzfile,'LEVEL',2);
ms2_scan_data = mzXML_struct_ms2.scan;
% ====> clear mzXML_struct_ms2


% Define the fields to keep
fieldsToKeep = {'num', 'msLevel', 'peaksCount', 'retentionTime', ...
               'startMz', 'endMz', 'lowMz', 'highMz', ...
               'basePeakMz', 'basePeakIntensity', 'totIonCurrent', 'precursorMz', 'peaks', 'msInstrumentID'};

% Remove all other fields from the struct array
ms2_scan_data = rmfield(ms2_scan_data, setdiff(fieldnames(ms2_scan_data), fieldsToKeep));
ms2_scan_xlsx = ms2_scan_data;
scan_num2 = cell2mat({ms2_scan_data.num})';

% ###################### MS2 peak data + retention time in seconds ###################### %
% ###################### MS2 peak data + retention time in seconds ###################### %
% ###################### MS2 peak data + retention time in seconds ###################### %
%%%%% TABELIZE PEAKS DATA  into [INTENSITY] absorbance [as a function of] m/z [MASS]
[ms2_peaks,ms2_ret_time] = mzxml2peaks(mzXML_struct_ms2, 'Levels', 2);

for i = 1:length(ms2_scan_data)
    %% DATA needed to create readable  .xlsx  files
    ms2_scan_xlsx(i).retentionTime = ms2_ret_time(i);
    ms2_scan_xlsx(i).precursorMz = ms2_scan_xlsx(i).precursorMz.value;

    %% DATA needed to create readable  .mat  files
    ms2_scan_data(i).retentionTime = ms2_ret_time(i); 
    ms2_scan_data(i).peaks = ms2_peaks(i);
    ms2_scan_data(i).precursorMz = ms2_scan_data(i).precursorMz.value;
end


% ###################### DETERMINE NUMBER of m1 & m2 EVENTS ###################### %
% ###################### DETERMINE NUMBER of m1 & m2 EVENTS ###################### %
% ###################### DETERMINE NUMBER of m1 & m2 EVENTS ###################### %
% ###################### DETERMINE NUMBER of m1 & m2 EVENTS ###################### %
% ###################### DETERMINE NUMBER of m1 & m2 EVENTS ###################### %

sel_ms1 = scan_num1;
sel_ms2 = scan_num2;

event_count = sel_ms1(2)-sel_ms1(1);  %%% M1 + MS2 events
total_number_ms2_scans  = 1 : 1 : length(ms2_scan_data);
event_indices  = reshape(total_number_ms2_scans, event_count-1, length(ms1_scan_data))';

%%% =====> SAVE ALL ms1 and ms2 EVENTS in seperate  .MAT files
% if ~exist('data', 'dir')
%         mkdir('ms1_data');
%         mkdir('ms2_data');
% end
mkdir('ms1_data');
mkdir('ms2_data');



save("./ms1_data/ms1_scan_data.mat","ms1_scan_data");
save("./ms1_data/ms1_mzXML.mat","mzXML_struct_ms1");


save("./ms2_data/ms2_scan_data.mat","ms2_scan_data");
save("./ms2_data/ms2_mzXML.mat","mzXML_struct_ms2");



% ##################### MS1 data event scans ######################## % 
% ##################### MS1 data event scans ######################## % 
cd('D:\RAW\EXAMPLE\20241022_TEST\ms1_data')
temp = struct2table(ms1_scan_xlsx);
filenam = 'ms1_scan_table.xlsx';
writetable(temp,filenam,'Sheet',1)

% ##################### MS2 data event scans ######################## % 
% ##################### MS2 data event scans ######################## % 
cd('D:\RAW\EXAMPLE\20241022_TEST\ms2_data')

for i=1:1:event_count-1  %% MS2 event_counts
    varnam = ['ms2_', num2str(i), '_scan']
    eval([varnam ' = ms2_scan_data(event_indices(:,i));'])

    % =====> SAVE AS .MAT file
    save(varnam,varnam)
    
    % % % ====> CREATE table and save as .xlsx file
    varnam = ['ms2_', num2str(i), '_scan_table']
    eval([varnam ' = ms2_scan_xlsx(event_indices(:,i));'])
    filenam = [varnam ,  '.xlsx'];
    temp = eval([ 'struct2table(' varnam ');']);
    writetable(temp,filenam,'Sheet',1);
end

% ##################### ############## ######################## % 

    
toc



%%
% ################################################################################## %
% ################################################################################## %
% ################################################################################## %
% ################################################################################## %
% ################################################################################## %
% ################################################################################## %
% ################################################################################## %
% ################################################################################## %
% ################################################################################## %
% ################################################################################## %
% ################################################################################## %

%%
% ##################### PEAK DATA STORAGE ####################### %
% ##################### PEAK DATA STORAGE ####################### %
% ##################### PEAK DATA STORAGE ####################### %
% ##################### PEAK DATA STORAGE ####################### %
% ##################### PEAK DATA STORAGE ####################### %
% ##################### PEAK DATA STORAGE ####################### %

clc
cd('D:\RAW\EXAMPLE\20241022_TEST\ms1_data')
load('ms1_scan_data.mat','ms1_scan_data' )

format bank

% Determine the desired number of filename characters (e.g., 15 characters)
num_characters = 15;
total_number_ms1_scans = numel(ms1_scan_data);

% Initialize an empty cell array to store the filenames
filenames = cell(total_number_ms1_scans, 1); 

% Loop to create filenames
for i = 1:1:total_number_ms1_scans

    % Generate the base filename
    base_filename = ['ms1_peak_scan_', num2str(i, '%03d'), '.xlsx']; 

    % Pad the filename with spaces to reach the desired length
    filename_padded = [base_filename, blanks(num_characters - length(base_filename))]; 
   
    % Store the padded filename in the cell array
    filenames{i} = filename_padded;

    % select + store excel data: [INTENSITY] absorbance [as a function of] m/z [MASS]
    temp = ms1_scan_data(i).peaks(1,1);
    xlsx_file = char(filenames(i));
    writematrix(temp{1}, xlsx_file ,'Sheet',1)
end

% Display the generated filenames
disp(filenames)



%%
% ################################################################################## %
% ################################################################################## %
% ################################################################################## %
% ################################################################################## %
% ################################################################################## %
% ################################################################################## %
% ################################################################################## %
% ################################################################################## %
% ################################################################################## %
% ################################################################################## %
% ################################################################################## %


%%
% ###### Visualizing and Preprocessing Hyphenated Mass Spectrometry ######
% ###### Metabolite and Protein/Peptide Profiling ###### %

tic

cd('D:\RAW\EXAMPLE\20241022_TEST')
info = mzxmlinfo('Lemonfiness_s3_onbespoten.mzXML','NUMOFLEVELS',true)

% The |MZXMLREAD| function reads the XML document into a MATLAB structure.
% The fields |scan| and |index| are placed at the first level of the output
% structure for improved access to the spectral data. The remainder of the
% mzXML document tree is parsed according to the schema specifications.
mzXML_struct = mzxmlread('Lemonfiness_s3_onbespoten.mzXML','LEVEL',1)

% To facilitate the handling of the data, the |MZXML2PEAKS| function
% extracts the list of peaks from each scan into a cell array (|peaks]|)
% and their respective retention time into a column vector (|time|). You
% can extract the spectra of certain level by setting the |LEVEL| input
% parameter.
 
[peaks,time] = mzxml2peaks(mzXML_struct);
numScans = numel(peaks);
basePeakInt = [mzXML_struct.scan.basePeakIntensity]';



toc

%%
% ################################################################################## %
% ################################################################################## %
% ################################################################################## %
% ################################################################################## %
% ################################################################################## %
% ################################################################################## %
% ################################################################################## %
% ################################################################################## %
% ################################################################################## %
% ################################################################################## %
% ################################################################################## %

%%
k = numScans;
j = 1
v = VideoWriter('Lemonfiness_ms1_onbespoten.avi');
open(v);

    fig = figure(1);
    temp = peaks(j,1);
    temp = temp{1};
    x = temp(:,1);
    y = temp(:,2); 
    h = text(1550,1100000,num2str(j))
    p = plot(x,y);
    axis([0 1550 0 1200000]);

% update data during animation
for j=1:1:k
    temp = peaks(j,1);
    temp = temp{1};
    p.XData = temp(:,1);      % update X and Y data properties of line object
    p.YData = temp(:,2);
    h = text(1550,1100000,num2str(j))
    frame = getframe(gcf);
    writeVideo(v,frame);

  

    pause(0.01);
    delete(h)
end
close(v)

% ==================>  findobj(fig,'Type','text')
    

%%
% Loop to create filenames
for i = 1:1:10%numScans

    % select + store excel data: [INTENSITY] absorbance [as a function of] m/z [MASS]
    temp = peaks(i,1);
    temp = temp{1};
    x = temp(:,1);
    y = temp(:,2);
    plot(x,y)
    pause(0.5)

end

%%

max(basePeakInt)


peaks_fil = cell(numScans,1);
for i = 1:numScans
    h = peaks{i}(:,2) > (basePeakInt(i).*0.75);
    peaks_fil{i} = peaks{i}(h,:);
end