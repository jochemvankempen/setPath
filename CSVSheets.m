classdef CSVSheets
    
    properties
        
        bookname = '';
        sheetname = '';
        filepath = '';
        header = {};
        cells = {};
        precision = [];
        
        BookSheetDelimiter = '.';% for constructing Book/Sheet filenames
        VariableDelimiter = '_';% for constructing Book/Sheet names in workspace
        
        CellSeparator = ',';% delimiter between cells
        TextDelimiter = '"';% indicating strings
        FileExtension = 'csv';
                
    end
    
    methods
        
        function SS = CSVSheets(SS1,varargin)
            % CSVSheets constructor
            %
            % copy structure
            % SS = CSVSheets(CSVSheets,SS2)
            %
            % set fields
            % SS = CSVSheets(CSVSheets,'field1',[value1] ... )
            
            if nargin>1
                f1 = fieldnames(SS1);
                SS = SS1;
                if length(varargin)==1 && isstruct(varargin{1})
                    % copy fields of structure
                    DF2 = varargin{1};
                    f2 = fieldnames(DF2);
                    for i=1:length(f2)
                        if ~any(strcmp(f2{i},f1));continue;end
                        SS.(f2{i}) = DF2.(f2{i});
                    end
                elseif  rem(length(varargin),2)==0
                    % set fields
                    f2 = varargin(1:2:end);
                    v2 = varargin(2:2:end);
                    for j = 1:numel(SS1)
                        for i=1:length(f2)
                            if ~any(strcmp(f2{i},f1));continue;end
                            SS(j).(f2{i}) = v2{i};
                        end
                    end
                end
            end
            
        end
        
        function [nRows,nCols] = size(SS)
            % returns size of the cells field
            % [nRows,nCols] = size(SS)
            [nRows,nCols] = size(SS.cells);
        end
                      
        function SS = read(SS,filepath,isHeadRow)
            % load a file into a SpreadSheet object
            % SS = read(CSVSheet,filepath,isHeadRow)
            if nargin<3;
                isHeadRow=true;
                if nargin<2
                    filepath = '';
                end;end
            if isempty(filepath) && isempty(SS.filepath)
                SS = uigetfile(SS);
            elseif ~isempty(filepath) && isempty(SS.filepath)
                SS.filepath = filepath;
            end
            if ~exist(SS.filepath,'file');error('File doesn''t exist: %s !',SS.filepath);end
            
            fid = fopen(SS.filepath,'r');
            if fid==-1;error('Error opening file %s !',SS.filepath);end
            fprintf('reading %s ... ',SS.filepath);
            EndOfFile = false;
            RowCnt = 0;
            while ~EndOfFile
                cRow = fgetl(fid);
                if cRow==-1
                    EndOfFile = true;
                else
                    RowCnt = RowCnt+1;
                    SS.cells{RowCnt,1} = cRow;
                end
            end
            SS = charline2col(SS);
            if isHeadRow
                SS.header = SS.cells(1,:);
%                 cellfun('isempty',SS.header):
                SS.cells(1,:) = [];
            end
            fprintf('detect numerical ... ');
            if size(SS)>0
                SS = getNumbers(SS);
            end
            fclose(fid);
            fprintf('done\n');
        end

%         function SS = read(SS,filepath)
%             % load a file into a SpreadSheet object
%             % SS = read(SS,filepath)
%             if nargin<2
%                 filepath = '';
%             end
%             if isempty(filepath) && isempty(SS.filepath)
%                 SS = uigetfile(SS);
%             elseif ~isempty(filepath) && isempty(SS.filepath)
%                 SS.filepath = filepath;
%             end
%             if ~exist(SS.filepath,'file');error('File doesn''t exist: %s !',SS.filepath);end
%             
%             [SS.cells,SS.header] = mexReadCSV(SS.filepath,SS.CellSeparator,SS.TextDelimiter);
%             SS.precision = cell(size(SS.cells));
%             
%             % convert char numbers to double
%             NumArray = str2double(SS.cells);
%             isNum = ~isnan(NumArray);
%             NumCols = all(isNum | cellfun('isempty',SS.cells),1);
%             
%             % get precision
%             SS.precision = cell(size(SS.cells));
%             SS.precision(isNum) = strfind(SS.cells(isNum),'.');
%             SS.precision(cellfun('isempty',SS.precision)) = {[0]};
%             SS.precision = cell2mat(SS.precision);
%             SS.precision(isNum&SS.precision>0) = cellfun('length',SS.cells(isNum&SS.precision>0))-SS.precision(isNum&SS.precision>0);
%             SS.precision(:,~NumCols) = NaN;
%             
%             % put numbers in cells
%             for iCol = find(NumCols)
%                 SS.cells(:,iCol) = num2cell(NumArray(:,iCol));
%             end
%             
%         end
        
        function SS = load(SS,BookDir,BookName,SheetNames)
            % load several files as a "workbook"
            % SS = load(SS,BookDir,BookName,SheetNames)
            % BookDir ...... [char] directory
            % Bookname ..... [char] Bookname
            % Sheetnames ... [cell]
            
            if nargin==4
                % if there any input arguments given, replace information
                % in the SpreadSheet
                SS = getWorkbook(SS,BookDir,BookName,SheetNames);
            elseif nargin==1
                % loop thru SS and use filepath in object
                SS = file2Bookname(SS);
            end
            
            % load files
            for iSS = 1:length(SS)
                SS(iSS) = read(SS(iSS),SS(iSS).filepath);
            end
        end
        
        function SS = save(SS,filepath)
            error('under construction !!');
        end
        
        function SS = write(SS,filepath,FormatString,DoAppend)
            % write object to CSV file
            % SS = write(SS,filepath,FormatString,DoAppend)
            %
            % if SS.filepath==1 ...  write to prompt
            
            % check input
            if nargin<2 && isempty(SS.filepath)  
                SS.filepath = fullfile('.',sprintf('%s%s%s.%s',SS.bookname,SS.BookSheetDelimiter,SS.sheetname,'csv'));
            elseif nargin>=2 && ~isempty(filepath)
                SS.filepath = filepath;
            end
            if nargin<3 || isempty(FormatString)
                FormatString = getFormatString(SS);
            end
            if nargin<4
                DoAppend = false;
            end
            [nRows,nCols] = size(SS);
            
            % open the file to write to
            if SS.filepath==1
                fid = 1;
            else
                if DoAppend
                    fprintf('appending %s ... ',SS.filepath);
                    fid = fopen(SS.filepath,'a');
                else
                    fprintf('writing %s ... ',SS.filepath);
                    fid = fopen(SS.filepath,'w');
                end
            end
            
            % print header
            if ~DoAppend
                for iCol = 1:nCols-1
                    fprintf(fid,'%s%s',SS.header{iCol},SS.CellSeparator);
                end
                fprintf(fid,'%s\n',SS.header{end});
            end
            
            % print cells
            hwait = waitbar(0,'writing spreadsheet');
            for iRow = 1:nRows
                for iCol = 1:nCols
                    if isnan(SS.cells{iRow,iCol})%xxx
                        SS.cells{iRow,iCol} = [];%xxx
                    end
                    if ~isempty(SS.TextDelimiter) ...%xxx
                            && ischar(SS.cells{iRow,iCol}) ...
                            && ~isempty(strfind(SS.cells{iRow,iCol},SS.CellSeparator))
                        SS.cells{iRow,iCol} = sprintf('%s%s%s',SS.TextDelimiter,SS.cells{iRow,iCol},SS.TextDelimiter);
                    end
                    waitbar(((iRow-1)*nCols+iCol)/(nRows*nCols),hwait);
                end
                fprintf(fid,FormatString,SS.cells{iRow,:});%xxx
                fprintf(fid,'\n');
            end
            close(hwait);
            if fid~=1;fclose(fid);end
            fprintf('done\n');

        end
        
        function SS = setCell(SS,value,RowNr,ColLabel)
            if ischar(ColLabel)
                iCol = findcolumn(SS,ColLabel);
            else
                iCol = ColLabel;
            end
            SS.cells{RowNr,iCol} = value;
        end
        
        function FormatString = getFormatString(SS,varargin)
            [nRows,nCols] = size(SS);
            iNumCol = isNumericalColumn(SS);
            FormatString = '';
            for iCol = 1:nCols
                if iNumCol(iCol) &&  isempty(SS.precision)
                    FormatString = cat(2,FormatString,'%1.0f');
                elseif iNumCol(iCol) &&  ~isempty(SS.precision)
                    FormatString = cat(2,FormatString,['%1.' sprintf('%1.0f',max(SS.precision(:,iCol))) 'f']);
                else
                    FormatString = cat(2,FormatString,'%s');
                end
                if iCol~=nCols
                    FormatString = cat(2,FormatString,SS.CellSeparator);
                end
            end
        end
             
        function SS = print(SS,SheetName,FormatString)
            % print spreadsheet content to the screen
            if nargin>=2
                iSheet = findSheet(SS,SheetName);
            else
                iSheet = 1:length(SS);
            end
            for iSS = iSheet
                nCol = size(SS.cells,1);
                i = cellfun('isclass',SS.cells,'char');
                for iR = 1:nCol
                    fprintf(FormatString,SS.cells{iR,:});
                end
            end
        end
        
        function display(SS,iRow,HeaderLabel,FormatString)
            % print spreadsheet content to the screen
            if nargin==1
                disp(SS);
            else
                nCol = length(HeaderLabel);
                nRow = length(iRow);
                iCol = findcolumn(SS,HeaderLabel);
                for iR = 1:nRow
                    fprintf('RowNr. %d\n',iRow(iR));
                    for iC = 1:nCol
                        fprintf('%20s : ',HeaderLabel{iC});
                        fprintf(FormatString{iC},SS.cells{iRow(iR),iCol(iC)});
                        fprintf('\n');
                    end
                    fprintf('\n');
                end
            end
        end

        function [i,SS] = rowfilter(SS,varargin)
            % filter rows of cells
            % [i,SS] = rowfilter(SS,'HEADER1',[VALUE1],...)
            % If the value field is left empty, all non empty cells are
            % filtered
            i = true(size(SS.cells,1),1);
            iNumCol = isNumericalColumn(SS);
            if nargin==1;return;end
            for iFilt = 1:2:length(varargin)
                iCol = strcmp(varargin{iFilt},SS.header);
                if iNumCol(iCol)
                    if ~isempty(varargin{iFilt+1})% look for filter parameter
                        i = i & ismember(cat(1,SS.cells{:,iCol}),varargin{iFilt+1});
                    else % look for non empty cells
                        i = i & ~isnan(cat(1,SS.cells{:,iCol}));
                    end
                else
                    if ~isempty(varargin{iFilt+1})% look for filter parameter
                        i = i & ismember(SS.cells(:,iCol),varargin{iFilt+1});
                    else
                        i = i & ~cellfun('isempty',SS.cells(:,iCol));
                    end
                end
            end
            SS.cells = SS.cells(i,:);
        end
        
        function [i,SS] = rowfilterOR(SS,varargin)
            % filter rows of cells
            % [i,SS] = rowfilter(SS,'HEADER1',[VALUE1],...)
            % If the value field is left empty, all non empty cells are
            % filtered
            i = false(size(SS.cells,1),1);
            iNumCol = isNumericalColumn(SS);
            if nargin==1;return;end
            for iFilt = 1:2:length(varargin)
                iCol = strcmp(varargin{iFilt},SS.header);
                if iNumCol(iCol)
                    if ~isempty(varargin{iFilt+1})% look for filter parameter
                        i = i | ismember(cat(1,SS.cells{:,iCol}),varargin{iFilt+1});
                    else % look for non empty cells
                        i = i | ~isnan(cat(1,SS.cells{:,iCol}));
                    end
                else
                    if ~isempty(varargin{iFilt+1})% look for filter parameter
                        i = i | ismember(SS.cells(:,iCol),varargin{iFilt+1});
                    else
                        i = i | ~cellfun('isempty',SS.cells(:,iCol));
                    end
                end
            end
            SS.cells = SS.cells(i,:);
        end
        
        function [i,SS] = inRange(SS,Hd,Bds,includeEdges)
            % check if values are within bounds
            % [i,SS] = inRange(SS,Hd,Bds,includeEdges)
            i = true(size(SS.cells,1),1);
            iNumCol = isNumericalColumn(SS);
            if nargin<4;includeEdges=true;end
            iCol = strcmp(Hd,SS.header);
            cVal = cat(1,SS.cells{:,iCol});
            if includeEdges
                i = cVal>=Bds(1)&cVal<=Bds(2);
            else
                i = cVal>Bds(1)&cVal<Bds(2);
            end
        end
        
        function [i,SS] = outRange(SS,Hd,Bds,includeEdges)
            % check if values are out of bounds
            % [i,SS] = outRange(SS,Hd,Bds,includeEdges)
            i = true(size(SS.cells,1),1);
            iNumCol = isNumericalColumn(SS);
            if nargin<4;includeEdges=true;end
            iCol = strcmp(Hd,SS.header);
            cVal = cat(1,SS.cells{:,iCol});
            if includeEdges
                i = cVal<=Bds(1)|cVal>=Bds(2);
            else
                i = cVal<Bds(1)|cVal>Bds(2);
            end
        end
        
        function VarName = checkout(SS)
            % assign to caller workspace for editing in variable editor
            % VarName = checkout(SS)
            for iSS = 1:length(SS)
                VarName{iSS} = sprintf('%s%s%s',SS(iSS).bookname,SS(iSS).VariableDelimiter,SS(iSS).sheetname);
                assignin('caller',VarName{iSS},cat(1,SS(iSS).header,SS(iSS).cells));
            end
        end
        
        function SS = checkin(SS,X,SheetName)
            % assign variables back to object after editing
            % SS = checkin(SS,X,SheetName)
            if nargin==3
                i = findSheet(SS,SheetName);
            else
                i =1;
            end
            SS(i).header = X(1,:);
            SS(i).cells = X(2:end,:);
        end
        
        function SS = charline2col(SS)
            % separate the lines of strings into cells
            % considers >CellSeparator< and >SS.TextDelimiter<
            ColAllocN = 300;
            [nRows] = size(SS.cells,1);
            X = cell(nRows,ColAllocN);
            isString = false(nRows,ColAllocN);
            nCols = zeros(nRows,1);
            for iRow=1:nRows
                % find delimiter that separates cells
                iSep = strfind(SS.cells{iRow},SS.CellSeparator);
                SepLength = length(SS.CellSeparator);
                nSep = length(iSep);
                % find delimiter that indicates text
                nTxtDel = 0;
                if ~isempty(SS.TextDelimiter)
                    TxtDelLength = length(SS.TextDelimiter);
                    iTxtDel = strfind(SS.cells{iRow},SS.TextDelimiter);
                    nTxtDel = length(iTxtDel);
                    if nTxtDel>0 && rem(nTxtDel,2)~=0
                        error('Detected and odd number of text/string delimiter!');
                    elseif nTxtDel>0 && nSep>0
                        iSepNot = false(1,nSep);
                        for j=1:nTxtDel/2
                            iSepNot(iSep>iTxtDel(j*2-1) & iSep<iTxtDel(j*2)) = true;
                        end
                        iSep(iSepNot) = [];
                    else
                        iSep = iSep;
                    end
                end
                nCols(iRow) = length(iSep)+1;
                
                % extract strings
                iSep = cat(2,1-SepLength,iSep,length(SS.cells{iRow})+1);
                for j=1:length(iSep)-1
                    X{iRow,j} = SS.cells{iRow}(iSep(j)+SepLength:iSep(j+1)-1);
                    % check if there are any Text delimiters in this cell
                    if nTxtDel>0&&any(iTxtDel>iSep(j)&iTxtDel<iSep(j+1))
                        cTxDel = strfind(X{iRow,j},SS.TextDelimiter);
                        X{iRow,j} = X{iRow,j}(cTxDel(1)+TxtDelLength:cTxDel(2)-1);
                    end
                end
            end
            X = X(:,1:max(nCols));
            SS.cells = X;
        end
        
        function SS = getNumbers(SS)
            % convert char numbers to double
            NumArray = str2double(SS.cells);
            isNum = ~isnan(NumArray);
            NumCols = all(isNum | cellfun('isempty',SS.cells),1);
            
            % get precision
            SS.precision = cell(size(SS.cells));
            SS.precision(isNum) = strfind(SS.cells(isNum),'.');
            SS.precision(cellfun('isempty',SS.precision)) = {[0]};
            SS.precision = cell2mat(SS.precision);
            SS.precision(isNum&SS.precision>0) = cellfun('length',SS.cells(isNum&SS.precision>0))-SS.precision(isNum&SS.precision>0);
            SS.precision(:,~NumCols) = NaN;
            
            % put numbers in cells
            for iCol = find(NumCols)
                SS.cells(:,iCol) = num2cell(NumArray(:,iCol));
            end
        end
        
        function i = isNumericalColumn(SS)
            i = cellfun('isclass',SS.cells,'double');
            i = all(i,1);
        end
        
        function i = isCharColumn(SS)
            i = cellfun('isclass',SS.cells,'char');
            i = all(i,1);
        end
        
        function i = isemptyCells(SS)
            [nRows,nCols] = size(SS);
            i = false(nRows,nCol);
            iNumCol = isNumericalColumn(SS);
            for iCol = 1:nCols
                if iNumCol(iCol)
                    i(:,iCol) = isnan(cat(1,SS.cells{:,iCol}));
                else
                    i(:,iCol) = cellfun('isempty',SS.cells(:,iCol));
                end
            end
        end
        
        function SS = uigetfile(SS,cDir)
            if nargin<2|| isempty(cDir)
                cDir = cd;
            end
            [filename, pathname, filterindex] = uigetfile( ...
                fullfile(cDir,'*.*'), ...
                'Pick a file', ...
                'MultiSelect', 'on');
            if pathname==0
                SS.filepath = [];
            elseif ischar(filename)
                SS.filepath = fullfile(pathname,filename);
            elseif iscell(filename)
                for i=1:length(filename)
                    SS(i).filepath = fullfile(pathname,filename{i});
                end
            end
        end
        
        function SS = getWorkbook(SS,BookDir,BookName,SheetNames)
            % get files that constitutes a workbook
            % BookDir ...... [char] directory
            % Bookname ..... [char] Bookname
            % Sheetnames ... [cell]            
            if isempty(BookName)
                % manually load file names
                SS = uigetfile(SS,BookDir);
                SS = file2Bookname(SS);
            elseif ~isempty(BookName) && isempty(SheetNames)
                % get sheetnames from filenames found
                if isempty(BookDir)
                    % ask for directory
                    BookDir = uigetdir(cd);
                    if BookDir==0;return;end
                end
                % search for files matching the bookname
                if length(SS)>1;error('There are more than one SpreadSheet object!');end
                cDir = dir(fullfile(BookDir,sprintf('*.%s',SS.FileExtension)));
                nCurrFileType = length(cDir);
                fNames = cell(nCurrFileType,1);
                [fNames{:}] = deal(cDir(:).name);
                fIndex = strmatch(BookName, fNames);
                if isempty(fIndex)
                    error('Can''t find bookname >%s<', SS.bookname);
                end
                for i=1:length(fIndex)
                    SS(i).filepath = fullfile(BookDir,fNames{fIndex(i)});
                    [fDir,fName,fExt] = fileparts(fNames{fIndex(i)});
                    [SS(i).bookname,SS(i).sheetname] = strread(fName,'%s%s','delimiter',SS(i).BookSheetDelimiter);
                    SS(i).bookname = char(SS(i).bookname);
                    SS(i).sheetname = char(SS(i).sheetname);
                end                
            elseif ~isempty(BookName) && ~isempty(SheetNames)
                if isempty(BookDir)
                    % ask for directory
                    BookDir = uigetdir(cd);
                    if BookDir==0;return;end
                end
                for i=1:length(SheetNames)
                    SS(i).bookname = BookName;
                    SS(i).sheetname = SheetNames{i};
                    SS(i).filepath = fullfile(BookDir,sprintf('%s%s%s.%s',BookName,SS(i).BookSheetDelimiter,SheetNames{i},SS(i).FileExtension));
               end
            end
        end
        
        function i = findSheet(SS,SheetName)
            nSh = length(SS);
            ShNms = cell(1,nSh);
            [ShNms{:}] = deal(SS(:).sheetname);
            if nargin<2
                i = ShNms;
            else
               i = find(strcmp(SheetName,ShNms));
            end
        end
        
        function SS = file2Bookname(SS)
            for i=1:length(SS)
                [fDir,fName,fExt] = fileparts(SS(i).filepath);
                [SS(i).bookname,SS(i).sheetname] = strread(fName,'%s%s','delimiter',SS(i).BookSheetDelimiter);
                SS(i).bookname = char(SS(i).bookname);
                SS(i).sheetname = char(SS(i).sheetname);
            end
        end
        
        function s = convert2struct(SS,idx)
            if nargin<2 || isempty(idx)
                idx = true(size(SS.cells,1),1);
            end
            for i=1:length(SS.header)
                if all(cellfun('isclass',SS.cells(idx,i),'double'))
                    s.(SS.header{i}) = cat(1,SS.cells{idx,i});
                else
                    s.(SS.header{i}) = SS.cells(idx,i);
                end
            end
        end
        
        function SS = insertcells(SS,c,cHeader,StartCell,ColPrecision)
            % insert cells, set precision of numerical cells
            [iRows,iCols] = size(SS);
            iNumCols = isNumericalColumn(SS);
            iCharCols = isCharColumn(SS);
            
            [icRows,icCols] = size(c);
            iNumCols_c = cellfun('isclass',c,'double');
            iCharCols_c = cellfun('isclass',c,'char');
            if any(diff(iNumCols_c,1,1)~=0,1) | any(diff(iCharCols_c,1,1)~=0,1)
                error('Cells to append have inconsistent num/char format!!');
            else
                iNumCols_c = all(iNumCols_c,1);
                iCharCols_c = all(iCharCols_c,1);
            end
            
            % get column indices
            if nargin>2
                [isHead,iHead] = ismember(cHeader,SS.header);
            else
                iHead = 1:icCols;
            end
            
            SS.cells(iRows+1:iRows+icRows,iHead) = c;
            
            % get precision
%             if  ~all((iNumCols(iHead) & iNumCols_c) | (iCharCols(iHead) & iCharCols_c))
%                 error('Please match cells to existing formats!');
%             end

            if nargin==5 && ~isempty(ColPrecision)
                SS.precision(iRows+1:iRows+icRows,iHead(iNumCols_c)) = repmat(ColPrecision(iNumCols_c),icRows,1);
            else
                if iRows==0 && isempty(SS.precision)
                    SS.precision = ones(icRows,iCols).*NaN;
                    SS.precision(:,iNumCols_c) = 0;
                elseif iRows>0 && ~isempty(SS.precision)
                    SS.precision(iRows+1:iRows+icRows,:) = repmat(max(SS.precision,[],1),icRows,1);% default is max precision of existing columns
                else
                    error('Need precision ..!!');
                end 
            end
                
        end
        
        function c = getcells(SS,ColumnName,RowIndex,OutputMode)
            % retrieve elements of the sheet
            % c = getcells(SS,ColumnName,RowIndex,OutputMode)
            %
            % OutputMode ... 'EVAL'
            %                'MAT'
            %                'CHAR'
            
            [nRow,nCol] = size(SS.cells);
            [ColFound,ColIndex] = ismember(ColumnName,SS.header);
            if nargin<3 || (nargin==4 && isempty(RowIndex))
                RowIndex = 1:nRow;
            elseif isempty(RowIndex) || ~any(RowIndex)
                c = {};
                return;
            end
            
            c = cell(length(find(RowIndex)),length(ColFound));
            c(:,ColFound) = SS.cells(RowIndex,ColIndex(ColFound));
            
            if nargin==4
                switch OutputMode
                    case 'EVAL'
                        if numel(c)==1 && ~isempty(c{1}) && ischar(c{1})
                            c = eval(char(c));
                        elseif numel(c)>1
                            okcells = cellfun('isclass',c,'char') & ~cellfun('isempty',c);
                            c(~okcells) = {[]};
                            for i = find(okcells(:)')
                                c{i} = eval(c{i});
                            end
                        end
                    case 'MAT'
                        if numel(c)==1 && ~isempty(c{1}) && isnumeric(c{1})
                            c = cell2mat(c);
                        elseif numel(c)>1
                            okcells = cellfun('isclass',c,'double') & cellfun('length',c)==1;
                            c(~okcells) = {NaN};
                            c = cell2mat(c);
                        else
                            c = NaN;
                        end
                    case 'CHAR'
                        if numel(c)==1 && ~isempty(c{1}) && ischar(c{1})
                            c = char(c);
                        elseif numel(c)>1
                            okcells = cellfun('isclass',c,'char') & cellfun('length',c)>0;
                            c(~okcells) = {''};
                            c = char(c);
                        else
                            c = NaN;
                        end
                end
            end
            
        end
        
        function [RowIndex,GroupedIndex] = listdlg(SS,ColHead,FilterPar,varargin)
            % use the listdlg function to get indices
            % [RowIndex,GroupedIndex] = listdlg(SS,ColHead,FilterPar,varargin)
            %
            % ColHead ..... construct list string from these headers
            % FilterPar ... header/value pairs for filtered list
            % varargin are properties for listdlg function
            %
            % RowIndex ....... vector of all row-indices
            % GroupedIndex ... row-indices grouped according to items shown
            %                  in the listdlg window
            RowIndex = [];
            GroupedIndex = {};
            
            n = size(SS.cells,1);
            if nargin<3 || isempty(FilterPar)
                FilterIndex = [1:n]';
            elseif ~isempty(FilterPar) && isnumeric(FilterPar)
                FilterIndex = FilterPar;
            else
                FilterIndex = find(rowfilter(SS,FilterPar{:}));
            end
            [ColFound,ColIndex] = ismember(ColHead,SS.header);
            nCol = length(ColFound);
            CellList = SS.cells(FilterIndex,ColIndex);
            
            % create formatstring
            NumCol = find(isNumericalColumn(SS));
            ListFormatString = '';
            for i=1:nCol
                if i>1
                    ListFormatString = [ListFormatString ' - '];
                end
                if any(NumCol==ColIndex(i))
                    ListFormatString = [ListFormatString '%1.0f'];
                else
                    ListFormatString = [ListFormatString '%s'];
                end
            end
            
            % create string for use with listdlg
            n = size(CellList,1);
            for i = 1:n;
                ListString{i,1} = sprintf(ListFormatString,CellList{i,:});
            end
            
            % check for duplicates within list, i.e. rows that have same
            % entries according to ColHead
            [uListString,uRowIndex,uListIndex] = unique(ListString);
            
            % apply listdlg
            [ListIndex,OK] = listdlg('ListString',uListString,varargin{:});
            if OK==0
                return;
            else
                RowIndex = ismember(uListIndex,ListIndex);% convert to filter index
                RowIndex = FilterIndex(RowIndex);% convert to row index
                if nargout>1
                    nList = length(ListIndex);
                    GroupedIndex = cell(nList,1);
                    for i=1:nList
                        GroupedIndex{i} = FilterIndex(uListIndex==ListIndex(i));
                    end
                end
            end
        end
        
        function [index,ColIndex] = find(SS,ColHead,ColElements,CaseInsensitive)
            % find rows that matches all the properties given in ColHead
            % index = find(SS,ColHead,ColElements,CaseInsensitive)
            % ColHead ........ {'Head1' 'Head2'}
            % ColElements .... [{Arr1} {Arr2}] Arr is cell or num, length
            %                  Arr indicates number of searches
            [ColFound,ColIndex] = ismember(ColHead,SS.header);
            nSearch = unique(cellfun('length',ColElements));
            if length(nSearch)>1
                error('Must be same number elements for each column header!');
            end
            m = length(ColHead);
            i = false(size(SS.cells,1),m,nSearch);
            for j=1:m
                if iscell(ColElements{j})
                    for k=1:nSearch
                        if CaseInsensitive
                            i(:,j,k) = strcmpi(ColElements{j}{k},SS.cells(:,ColIndex(j)));
                        else
                            i(:,j,k) = strcmp(ColElements{j}{k},SS.cells(:,ColIndex(j)));
                        end
                    end
                elseif isnumeric(ColElements{j})
                    for k=1:nSearch
                        i(:,j,k) = ColElements{j}(k)==cat(1,SS.cells{:,ColIndex(j)});
                    end
                else
                    error('Search terms should be either a cell array or numerical!');
                end
            end
%             if any(sum(all(i,2))>1)
%                 error('Found more than one index per element!');
%             end
            [index,dummy] = find(all(i,2));
        end
                
        function [B,I,J,G] = unique(SS,ColHead,RowIndex)
            % Finds unique rows within columns defined by ColHead.
            % B ... are the unique cells
            % I ... A(I) = B
            % J ... B(J) = A
            % RowIndex ... apply only on these rows
            % G ... each cell contains row indices for one element in B
            [ColFound,ColIndex] = ismember(ColHead,SS.header);
            NumCol = find(isNumericalColumn(SS));
            nCol = length(ColHead);
            if nargin<3
                RowIndex = 1:size(SS.cells,1);
            elseif isempty(RowIndex)
                B = [];I = [];J = [];G = [];
                return;
            end
            nRows = length(RowIndex);
            B = SS.cells(RowIndex,ColIndex);
            I = ones(nRows,1);
            J = ones(nRows,1);
            RedundantRow = false(nRows,1);
            % go thru each row and see if you find the same row in the rows
            % before, if you do clear row content
            for i = 2:nRows
                isNotCell = cellfun('isempty',B(1:i-1,:));
                isSame = false(i-1,nCol);
                for iCol = 1:nCol
                    if any(NumCol==ColIndex(iCol))
                        isSame(~isNotCell(:,iCol),iCol) = cat(1,B{1:i-1,iCol}) == B{i,iCol};
                    else
                        isSame(:,iCol) = strcmp(B{i,iCol},B(1:i-1,iCol));
                    end
                end
                if any(all(isSame,2))% find same rows
                    B(i,:) = {[]};
                    I(i,1) = find(all(isSame,2));% which one is the same row
                    RedundantRow(i,1) = true;
                    J(i,1) = sum(~RedundantRow(1:I(i,1),1));
                else
                    I(i,1) = i; 
                    J(i,1) = sum(~RedundantRow(1:i,1));
                end
            end
            B = B(~RedundantRow,:);
            I = I(~RedundantRow);
            I = RowIndex(I);%
            % get grouped indices
            for i = 1:size(B,1)
                G{i,1} = RowIndex(J==i);
            end
        end
        
        function i = findcolumn(SS,ColLabel)
            if ischar(ColLabel)
                i = find(strcmp(ColLabel,SS.header));
            elseif iscell(ColLabel)
                [isFound,index] = ismember(ColLabel,SS.header);
                i = zeros(size(ColLabel));
                i(isFound) = index(isFound);
                i(~isFound) = NaN;
            end
        end
        
        function SS = extractRows(SS,i)
            SS.cells = SS.cells(i,:);
        end
        
        function SS = datenum(SS,DateColumn,varargin)
            if isnumeric(DateColumn)
                dci = DateColumn;
            else
                dci = findcolumn(SS,DateColumn);
            end
            x = datenum(SS.cells(:,dci),varargin{:});
            SS.cells(:,dci) = num2cell(x);
            
        end
            
        function SS = datestr(SS,DateColumn,varargin)
            if isnumeric(DateColumn)
                dci = DateColumn;
            else
                dci = findcolumn(SS,DateColumn);
            end
            x = datestr(cat(1,SS.cells{:,dci}),varargin{:});
            SS.cells(:,dci) = cellstr(x);
            
        end
            
        function SS = insertRows(SS,i,indexOffset,n)
            [nRows,nCols] = size(SS.cells);
            if isnumeric(indexOffset)
                i = i+1;
            elseif ischar(indexOffset) && strcmpi(indexOffset,'below')
                i = i+1;
            elseif ischar(indexOffset) && strcmpi(indexOffset,'above')
                i = i-1;
            end
            SS.cells = cat(1, ...
                SS.cells(1:i-1,:), ...
                cell(n,nCols), ...
                SS.cells(i:end,:));
        end
        
        function SS = copyRow(SS,sourceRow,targetRow,ColumnHeaders)
            iCols = findcolumn(SS,ColumnHeaders);
            SS.cells(targetRow,iCols) = SS.cells(sourceRow,iCols);
        end
        
    end
end