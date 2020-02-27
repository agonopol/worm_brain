function [indentation] = indent(indentationLevel)
    indentation = '';
    for i=1:indentationLevel
        indentation = [indentation '  '];
    end
end
