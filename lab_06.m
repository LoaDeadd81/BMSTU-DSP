function lab_06
    mult = 5;
    step = 0.005;
    t = -mult:step:mult;
    
    A = 1.0;
    sigma = 0.5;
    x0 = gauspls(t, A, sigma);
    
    NA = 0;
    NS = 0.05;
    n1 = normrnd(NA, NS, [1 length(x0)]);
    x1 = x0 + n1;
    
    count = 7;
    M = 0.4;
    n2 = impnoise(length(x0), count, M);
    x2 = x0 + n2;
    
    y1 = wiener(fft(x1), fft(n1));
    y2 = wiener(fft(x2), fft(n2));

    figure(1)
    plot(t,x1);
    title('Помеха по Гауссу');

    figure(2)
    plot(t,x2);
    title('Помеха по равномерному закону');
    
    figure(3)
    plot(t, ifft(fft(x1).*y1));
    title('Помеха по Гауссу с фильтром Винера');
    
    figure(4)
    plot(t, ifft(fft(x2).*y2));
    title('Помеха по равномерному закону с фильтром Винера');

    function y = gauspls(x,A,s)
	    y = A * exp(-(x / s).^2);
    end

    function y = impnoise(size, N, mult)
        step = floor(size / N);
        y = zeros(1, size);
        for i = 1:floor(N / 2)
            y(round(size / 2) + i * step) = mult * (0.5 + rand);
            y(round(size / 2) - i * step) = mult * (0.5 + rand);
        end
    end

    function y = wiener(x,n)
        y = 1 - (n./x).^2;
    end
end