function heatmaptopannot(annot, yoffset, fontsize)

    n = length(annot);
    for i = 1:n
        h = text(i, yoffset, annot{i}, 'fontsize', fontsize);
        set(h, 'rotation', 90);
    end