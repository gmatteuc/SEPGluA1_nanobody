function [sortedNames, ndx] = sort_nat(cellstr_in)
% sort_nat  Natural‐order sort of filenames containing "chanNN_"
%   [sortedNames,ndx] = sort_nat(cellstr_in)
%   Extract the NN after "chan", sort numerically, and return
%   the sorted list plus the index vector.

    nums = Inf(size(cellstr_in));
    for i = 1:numel(cellstr_in)
        tok = regexp(cellstr_in{i}, 'chan(\d+)_', 'tokens');
        if ~isempty(tok)
            nums(i) = str2double(tok{1}{1});
        end
    end
    [~, ndx]    = sort(nums);
    sortedNames = cellstr_in(ndx);
end

