%% Colorimetrical data (CIEXYZ or CIELAB) of ColorChecker
% Measurement data of the 24/30/48 color patches are stored as a class.
% 12/1/2023

% called by cprapp with LAB_Reference and LAB_Device:
%   setLAB()

classdef DeviceData < handle
    
    properties
        source = ''                                                         % record the source data color space: 'spd', 'XYZ', or 'LAB'
         
        Lab = zeros(24,3)                                                   % CIELAB data
        spec                                                                % spectral data
        XYZ                                                                 % CIEXYZ of the 24 patches -- 24x3
        XYZ_n                                                               % CIEXYZ of the white patch -- 1x3
        rgb                                                                 % sRGB data -- 24x3
    end
    
    methods

        function setXYZ (obj, XYZ)
            % Enter the data measured in CIEXYZ
            % called by ColorPerformanceReview::one_pager_from_XYZ()            

            % XYZ is 25x3 for 24-patch ColorChecker
            %assert(size(XYZ,1)==25)                                        % to-do: check the data size

            n = size(XYZ,1);

            obj.XYZ = XYZ(1:n-1,:);
            obj.XYZ_n = XYZ(n,:);

            for i = 1:n-1
                obj.Lab(i,:) = CPR.ColorClass.XYZ2lab(obj.XYZ(i,:),obj.XYZ_n);
            end

            for i = 1:n-1
                rgb = lab2rgb(obj.Lab(i,:));
                rgb = min(1,rgb);
                rgb = max(0,rgb);
                obj.rgb(i,:) = uint8(rgb*255);
            end

            obj.source = 'XYZ';
        end

       function setLAB (obj, LAB)
            % Enter the data measured in CIELAB
            % LAB is 24x3

            n = size(LAB,1);

            %assert(size(RGB,1)==24)

            obj.Lab = LAB(1:n,:);

            for i = 1:n
                obj.Lab(i,:) = LAB(i,:);
            end

            obj.source = 'LAB';
        end

         
        function setRGB (obj, RGB)
            % Enter the data measured in sRGB
            % RGB is 24x3

            n = size(RGB,1);

            %assert(size(RGB,1)==24)

            obj.rgb = RGB(1:n,:);

            for i = 1:n
                obj.Lab(i,:) = rgb2lab(obj.rgb(i,:)/255,'colorspace','srgb','whitepoint','d65');
            end

            obj.source = 'RGB';
        end

        function calcRGB (obj)
            % Convert the CIELAB to sRGB for visualization

            for i = 1:size(obj.Lab,1)
                rgb = lab2rgb(obj.Lab(i,:));

                rgb = min(1,rgb);                                           % clip the data
                rgb = max(0,rgb);

                obj.rgb(i,:) = uint8(rgb*255);
            end
            
        end
        
    end
end

