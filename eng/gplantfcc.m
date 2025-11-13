function fig = gplantfcc(varargin)
% gplantfcc(...)
% FCC
ffig = ne_group(varargin,'FCC','pplantfccf','pplantfccflow1','pplantfccflow2','pplantfcct','pplantfccs');
if nargout > 0 fig = ffig; end
