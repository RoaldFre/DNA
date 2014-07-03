% Load the native data from the file and save it again as float-binary with 
% max gzip compression
function toFloatBinary(file)

try
	data = load(file);
catch
	disp(["Couldn't load file ",file," BAILING OUT!"]);
	return
end

save("-float-binary", [file,"_floatBinary"], '-struct', 'data');

% manually compress with best compression instead of giving '-z' flag in 
% save() command
[ret,_] = system(['gzip --best ',file,'_floatBinary && mv ',file,'_floatBinary.gz ',file]);

if ret
	disp(["Couldn't zip or move file ",file," BAILING OUT!"]);
	return
end

disp(["Successfully converted file ",file," :-)"]);

