% Load the ascii matrix from the file, and save it to Octave's native 
% binary format, overwriting the original file. The matrix is stored in the 
% variable 'data' of the resulting file.
%
% All comments in the file are saved to the 'comments' variable. A comment 
% is a line that starts with a '#'.
function toNative(file)

[_, comments] = system(["grep '^#' ",file]);
if not(isempty(comments))
	% Cut off last newline
	comments = comments(1 : end-1);
end
data = load(file);
save("-binary", file, "data", "comments");
% manually compress with best compression instead of giving '-z' flag in 
% save() command
system(['gzip --best ',file,'; mv ',file,'.gz ',file]);

