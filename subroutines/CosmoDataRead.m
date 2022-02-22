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
        [num,~,~]    = xlsread(filename,find(nuclide_sheet));          % load CRN data
        [Xdata,Xtxt,Xraw] = xlsread(filename,find(comp_sheet));        % load compositional data,
        Xraw = Xraw(2:end,:);                                          % remove header     
        
        if and(isnan(Xraw{ind,2}),~isnan(Xraw{ind,5}))  % check if data are enterd correctly
            X.mode = 'soil';             
        elseif and(isnan(Xraw{ind,5}), ~isnan(Xraw{ind,2}))
            X.mode = 'bedrock'; 
        elseif X.n == 1 && and(isnan(Xraw{ind,5}), isnan(Xraw{ind,2}))
            X.mode = 'noComp'; 
        else
            error('Please check your input table. You did not provide the data as specified in the example input')
        end
            

        X.W              = Xraw{ind,11};                                    % Weathering rate
        X.Wstd           = Xraw{ind,12};
        X.soil_mass      = Xraw{ind,8};
        
        num = num(ind,:); 
        sampName = Xraw{ind,1}; 
        
    case 3          % paired CRN
        [num10,~,~] = xlsread(filename,'10Be Cronus');               % load 10Be data
        [num36,~,~]     = xlsread(filename,'36Cl Cronus');               % load 36Cl data
        [Xdata,Xtxt,Xraw]  = xlsread(filename,find(~nuclide_sheet));     % load compositional data
        Xraw = Xraw(2:end,:);                                            % remove header 
        
        if size(Xdata,2)== 4
            X.mode = 'soil';       
        elseif size(Xdata,2)== 7
            X.mode = 'bedrock';
        end
        
        num10 = num10(ind,:);   
        num36 = num36(ind,:);
        sampName = Xraw{ind,:}; 
        num.num10 = num10; num.num36 = num36;
        X.soil_mass    = Xraw{ind,8};
end

DEMdata.method = Xraw{ind,10};
global scaling_model
scaling_model  = Xraw{ind,9};       


switch X.mode
    case 'soil'
        X.fQzS = Xraw{ind,5}; 
        X.fCaS = Xraw{ind,6}; 
        X.fXS  = Xraw{ind,7};  
    case 'bedrock'
        X.fQzB = Xraw{ind,2}; 
        X.fCaB = Xraw{ind,3}; 
        X.fXB  = Xraw{ind,4};  
end


end

