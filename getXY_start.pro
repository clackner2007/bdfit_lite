;+
;gets the xy start point by getting the brightest pixels
;
;-
FUNCTION getXY_start, image
  xlen = (size(image, /dimensions))[0]
  ylen = (size(image, /dimensions))[1]
  ys = reverse(sort(smooth(total(image,1),min([4, ylen-1]))))
  xs=reverse(sort(smooth(total(image,2),min([4, xlen-1]))))
  y0=ys[0]
  x0=xs[0]
  return, [x0, y0]
END
