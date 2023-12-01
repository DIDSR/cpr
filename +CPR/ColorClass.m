%% Class for calculating color differences in CIE dE76, dE94, and CIEDE2000
% Each object is one color

classdef ColorClass

    %%
    properties

        lab = [0 0 0];            % store color in CIELAB

    end

    %%
    %  Methods for the objects
    %
    methods

        function obj = ColorClass (lab)
            % initialize the color with CIELAB

            obj.lab = lab;
        end

        function info (obj)
            % print the CIELAB values

            fprintf('L*a*b* = (%.2f,%.2f,%.2f)\n',obj.lab(1),obj.lab(2),obj.lab(3));
        end

        function ret = labget (obj)
            % get the CIELAB values as a 1x3 vector

            ret = [obj.lstar obj.astar obj.bstar];
        end

        function dE00 = minus (obj1, obj2)
            % override the '-' operator with dE

            dE00 = CPR.ColorClass.lab2dE00(obj1.lab,obj2.lab);
        end

    end

    %%
    %  Static methods for public access
    %
    methods (Static)

        function lab = XYZ2lab (XYZ, XYZ_white)
            % convert CIEXYZ to CIELAB w.r.t to the reference white

            lab = xyz2lab(XYZ,'WhitePoint',XYZ_white);
        end

        function [dE00 dE94 dEab] = LAB2dE (lab1, lab2)
            % calculate color differene between _lab1_ and _lab2_ in dE00, dE94, and dEab
          
            dE00 = CPR.ColorClass.lab2dE00(lab1', lab2');
            dE94 = CPR.ColorClass.lab2dE94(lab1', lab2');
            dEab = CPR.ColorClass.lab2dEab(lab1', lab2');
        end

        function dE00 = lab2dE00 (lab1, lab2)
            % calculate CIEDE2000

            dE00 = imcolordiff(lab1,lab2,'isInputLab',true,"Standard",'CIEDE2000');
        end

        function dE94 = lab2dE94 (lab1, lab2)
            % calcualte CIE dE 1994

            dE94 = imcolordiff(lab1,lab2,'isInputLab',true,"Standard",'CIE94');
        end
        
        function dEab = lab2dEab (lab1, lab2)
            % calculate CIE dE 1976

            dEab = deltaE(lab1,lab2,"isInputLab",true);
        end        

    end

end