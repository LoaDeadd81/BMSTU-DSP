function lab_05
    sigma = 0.5;

    mult = 5;
    step = 0.005;
    t = -mult:step:mult;

    x0 = gauspls(t, sigma);

    NM = 0;
    NS = 0.05;
    n1 = normrnd(NM,NS,[1 length(x0)]);
    x1 = x0+n1;

    count = 7;
    M = 0.4;
    n2 = impnoise(length(x0),count,M);
    x2 = x0+n2;

    G = gaussfilt(4,20);
    BB = buttfilt(6,20);

    figure(1)
    plot(t,x1);
    title('Помеха по Гауссу');

    figure(2)
    plot(t,x2);
    title('Помеха по равномерному закону');

    figure(3)
    plot(t,x1-filtfilt(G,1,x1));
    title('Помеха по Гауссу с фильтром Гауса');

    figure(4)
    plot(t,x1-filtfilt(BB,1,x1));
    title('Помеха по Гауссу с фильтром Баттеруорта');

    figure(5)
    plot(t,x2-filtfilt(G,1,x2));
    title('Помеха по равномерному закону с фильтром Гауса');

    figure(6)
    plot(t,x2-filtfilt(BB,1,x2));
    title('Помеха по равномерному закону с фильтром Баттеруорта');
end

function y = gauspls(x,s)
    y = exp(-(x/s).^2);
end

function y = impnoise(size,N,mult)
    step = floor(size/N);
    y = zeros(1,size);
    for i = 1:floor(N/2)
        y(round(size/2)+i*step) = mult*(0.5+rand);
        y(round(size/2)-i*step) = mult*(0.5+rand);
    end
end

function y = buttfilt(D,size)
    x = linspace(-size/2,size/2,size);
    y = 1./(1+(D./x).^4);
    y = y/sum(y);
end

function y = gaussfilt(sigma,size)
    x = linspace(-size/2,size/2,size);
    y = 1 - exp(-x.^2/(2*sigma^2));
    y = y/sum(y);
end