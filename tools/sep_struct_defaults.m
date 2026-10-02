function P = sep_struct_defaults(P, D)
%SEP_STRUCT_DEFAULTS  Fill the missing fields of an options struct from defaults.
%   P = SEP_STRUCT_DEFAULTS(P, D) copies every field of D that P does not
%   have. Fields already in P are left unchanged.
%
%   Example:
%       P = sep_struct_defaults(P, struct('max_detail', 5, 'verbose', true));

for f = fieldnames(D)'
    if ~isfield(P, f{1})
        P.(f{1}) = D.(f{1});
    end
end
end
