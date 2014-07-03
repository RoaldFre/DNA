% Load the ascii matrix from the file, and save it to Octave's native 
% binary format, overwriting the original file. The matrix is stored in the 
% variable 'data' of the resulting file.
%
% All comments in the file are saved to the 'comments' variable. A comment 
% is a line that starts with a '#'.
function toNative(file, singlePrecision)

if nargin < 2
	singlePrecision = false;
end

[_, comments] = system(["grep '^#' ",file]);
if not(isempty(comments))
	% Cut off last newline
	comments = comments(1 : end-1);
end
data = load(file);
if singlePrecision
	save("-float-binary", file, "data", "comments");
else
	save("-binary",       file, "data", "comments");
end
% manually compress with best compression instead of giving '-z' flag in 
% save() command
system(['gzip --best ',file,'&& mv ',file,'.gz ',file]);

