function outvar = scaleVar(invar, method, dim)
% scales matrix of values according to method
%
% jochem van kempen, 20-03-2018

% [x, y] = size(invar);

% if x == 1;
%     invar = invar';
% end
if nargin < 3
    dim = [];
end

switch method
    case 'minmax'
        %%% scales to the minimum and maximum value, so that it
        %%% will be scaled between 0 and 1
        %%% yi=(xi-min(xi)) /(max(xi)-min(xi))
%         outvar = (invar - repmat(min(invar(:)),length(invar),1)) ./ (repmat(max(invar(:)),length(invar),1) - repmat(min(invar(:)),length(invar),1));
        
        outvar = (invar - min(invar(:))) / (max(invar(:)) - min(invar(:)));

    case 'max'
        %%% scales by maximum value       
        outvar = invar / max(invar(:));
        
    case {'maxdim','minmaxdim'}
        %%% scales by maximum value within a given dimension   
        formula = [];
        for idim = 1:ndims(invar)
            if idim == dim
                formula = [formula 'i' ];
            else
                formula = [formula ':' ];
            end
            if idim ~= ndims(invar)
                formula = [formula ',' ];
            end
        end
        
        outvar = zeros(size(invar));
        for i = 1:size(invar,dim)
            
            eval(['tmp = squeeze(invar(' formula '));'])
            switch method
                case 'maxdim'
                    eval(['outvar(' formula ') = tmp / max(tmp(:));']);
                case 'minmaxdim'
                    eval(['outvar(' formula ') = (tmp / max(tmp(:))) / (max(tmp(:)) - min(tmp(:)));']);
            end
        end
            
end