function [num,sampName,X,DEMdata] = CosmoDataRead(filename,varargin)
% Read the information from the excel sheets. 
% If you have several samples in your input sheet provide the number of the 
% sample you want to run as second input argument.
% E.g. CosmoDataRead(filename,2)

[~,sheets] = xlsfinfo(filename);

if ~isempty(varargin)
    ind  = varargin{1}; 
else
    ind  = 1;
end

nuclide_sheet = ~cellfun(@isempty,strfind(sheets, 'Cronus'));    % find the excel sheet that has the nuclide information
comp_sheet = ~cellfun(@isempty,strfind(sheets, 'Composition'));  % find the excel sheet that has the compositional information

switch length(sheets)
    case 2           % single CRN plus Weathering
        nuclide = sheets{nuclide_sheet}(1:4);                            % identify which nuclide was provided
        switch nuclide, case '10Be', X.n = 1; case '36Cl', X.n = 2; end  % save which nuclide is being used  
        
        % load data
        [num,txt,~]    = xlsread(filename,find(nuclide_sheet));          % load CRN data
        [Xdata,Xtxt,~] = xlsread(filename,find(comp_sheet));         % load compositional data,
        
        if size(Xdata,2)== 8
            X.mode = 'soil';
        elseif size(Xdata,2)== 11
            X.mode = 'bedrock';
        else
            error('Please check your input table. You did not provide the data as specified in the example input')
        end
            
        switch X.mode
            case 'soil'
                X.W              = Xdata(ind,7);                                    % Weathering rate
                X.Wstd           = Xdata(ind,8);
                X.soil_mass      = Xdata(ind,4);
            case 'bedrock'
                X.W              = Xdata(ind,10);                                    % Weathering rate
                X.Wstd           = Xdata(ind,11);
                X.soil_mass      = Xdata(ind,7);
            end
        
        num = num(ind,:); 
        sampName = txt{ind+2,1}; 
        
    case 3          % paired CRN
        [num10,txt10,~] = xlsread(filename,'10Be Cronus');               % load 10Be data
        [num36,~,~]     = xlsread(filename,'36Cl Cronus');               % load 36Cl data
        [Xdata,Xtxt,~]  = xlsread(filename,find(~nuclide_sheet));        % load compositional data
        
        num10 = num10(ind,:);   
        num36 = num36(ind,:);
        sampName = txt10{ind+2,:}; 
        num.num10 = num10; num.num36 = num36;
        X.soil_mass    = Xdata(ind,7);
end

DEMdata.method = Xtxt{ind+1,10};
global scaling_model
scaling_model  = Xtxt{ind+1,9};       


switch X.mode
    case 'soil'
        X.fQzS = Xdata(ind,1); 
        X.fCaS = Xdata(ind,2); 
        X.fXS  = Xdata(ind,3);  
    case 'bedrock'
        X.fQzB = Xdata(ind,1); 
        X.fCaB = Xdata(ind,2); 
        X.fXB  = Xdata(ind,3);  
end


end

