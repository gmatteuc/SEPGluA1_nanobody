function h = sep_hash_array(x)
%SEP_HASH_ARRAY  MD5 hash of a numeric, logical or char array.
%   h = SEP_HASH_ARRAY(x) returns 'class:[size]:md5' as a char array. Two
%   arrays get the same hash only if class, size and every value are the same
%   (NaN included).

if isempty(x)
    h = sprintf('empty:%s:%s', class(x), mat2str(size(x)));
    return
end

% bytes of logical and char arrays via a numeric type
if islogical(x)
    x = uint8(x);
end
if ischar(x)
    x = uint16(x);
end
assert(isnumeric(x), 'sep_hash_array: unsupported class %s', class(x));
if ~isreal(x)
    x = [real(x(:)); imag(x(:))];
end

% feed the bytes in chunks (java arrays are limited in size); the digest is
% the same as for a single update
md = java.security.MessageDigest.getInstance('MD5');
n = numel(x);
chunk = 2^24;
for i0 = 1:chunk:n
    part = x(i0:min(n, i0 + chunk - 1));
    md.update(typecast(part(:), 'uint8'));
end
digest = typecast(md.digest(), 'uint8');
h = sprintf('%s:%s:%s', class(x), mat2str(size(x)), sprintf('%02x', digest));
end
