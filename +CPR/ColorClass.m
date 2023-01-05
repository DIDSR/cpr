%% Class for color
classdef ColorClass

    properties
        lab = [0 0 0];
    end

    methods

        function obj = ColorClass (lab)
            obj.lab = lab;
        end

        function info (obj)
            fprintf('L*a*b* = (%.4f,%.4f,%.4f)\n',obj.lstar,obj.astar,obj.bstar);
        end

        function ret = labget (obj)
            ret = [obj.lstar obj.astar obj.bstar];
        end

        function dE00 = minus (obj1, obj2)
            dE00 = ColorClass.lab2dE00(obj1.lab,obj2.lab);
        end

    end

    methods (Static)

        function lab = XYZ2lab (XYZ, XYZ_white)
            lab = xyz2lab(XYZ,'WhitePoint',XYZ_white);
        end

        function [dE00 dE94 dEab] = LAB2dE (lab1, lab2)
            dE00 = CPR.ColorClass.lab2dE00(lab1', lab2');
            dE94 = CPR.ColorClass.lab2dE94(lab1', lab2');
            dEab = CPR.ColorClass.lab2dEab(lab1', lab2');
        end

        function dE00 = lab2dE00 (lab1, lab2)
            dE00 = imcolordiff(lab1,lab2,'isInputLab',true,"Standard",'CIEDE2000');
        end

        function dE94 = lab2dE94 (lab1, lab2)
            dE94 = imcolordiff(lab1,lab2,'isInputLab',true,"Standard",'CIE94');
        end
        
        function dEab = lab2dEab (lab1, lab2)
            dEab = deltaE(lab1,lab2,"isInputLab",true);
        end        

    end

end