classdef ColorClassXYZ < handle

    properties
        XYZ
        XYZ_n
        xyz
        xyz_n
        lab = []
        rgb
        cct
    end

    methods

        function obj = ColorClassXYZ (XYZ,XYZ_n)

            assert(nargin >= 1)

            obj.XYZ = XYZ;
            obj.xyz = XYZ ./ sum(XYZ);
            obj.calculate_CCT;

            if nargin == 2
                obj.XYZ_n = XYZ_n;
                obj.xyz_n = XYZ_n ./ sum(XYZ_n);

                cc = ColorConversionClass;
                obj.lab = cc.XYZ2lab(obj.XYZ,obj.XYZ_n);
                obj.rgb = uint8(lab2rgb(obj.lab) * 255);
            end

        end

        function im = create_color_patch_xyz (obj, Y)
            xyz = obj.xyz;
            XYZ = xyz * Y / xyz(2);
            rgb = xyz2rgb(XYZ)*255;
            im = obj.create_color_patch_rgb(rgb);
        end

        function im = create_color_patch_rgb (obj,rgb)
            im = uint8(zeros(100,100,3));
            im(:,:,1) = rgb(1);
            im(:,:,2) = rgb(2);
            im(:,:,3) = rgb(3);
            %imshow(im);
        end

        function ret = calculate_CCT (obj)
            x = obj.xyz(1);
            y = obj.xyz(2);

            n = (x-0.3320)/(0.1858-y);
            ret = 437*n^3 + 3601*n^2 + 6861*n + 5517;
            obj.cct = ret;
        end
    end

end