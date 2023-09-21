%% Random Channels %%
% To generate random shaped channels, Equations are from Schmidt 2008
close all, clear all
reverse = false;
criteria = false;
while criteria == false 

    a = [2];% 30 20]; % Minimal width of the channel
    Dr = 0.5;% Dmin/Dmax 
    L0 = [15]; % Correlation length, irregularity of the channel % 10 for irregular
    L=150;
    z = linspace(0,L,300); % channel vertical location
    for i = 1:100
        e = rand;
        k = 1/2*(2*pi/L)*i;
        if k >= 2*pi/L0
%             display('Continue') 
            break
        end
        D(i,1:300) = (a/sqrt(1+(k*L0)^2) .* sin(k.* z + 2*pi.*e)); % Channel diameter
    end
    D = sum(D);
    D = D + abs(min(D));
    D = D + max(D);
    x = (Dr*max(D)-min(D))/(1-Dr);
    D = D + x;
    r = D/2;
    r_throat = min(r);
    throat = find(r_throat == r); %n-loaction of throat
    if length(throat)>1
        display("Warning: 2 locations for throat found")
    end
    throat = throat(1);
    z_throat = z(throat);
    figure(1);
    plot([z_throat z_throat],[r_throat -r_throat],'--')
    hold on
    plot(z,r,z,-r)
    hold off 
    ylabel('Radius')
    xlabel('z')
    d_e_d_s = r(end)/r_throat;
    if d_e_d_s > 1.42 && d_e_d_s < 1.437 && z_throat > 110.8 && z_throat < 111.2 % && r(1)> 1.2*r(end)
    %if d_e_d_s > 5 && d_e_d_s < 5.5 && z_throat > 110.8 && z_throat < 111.2 && r(1) < 0.8*r(end)    
        criteria = true
    end
end
%% Adjust final diameter %%
%r = [6.69901118328829,6.71029586260731,6.72008078968591,6.72833210261095,6.73501962811677,6.74011698621199,6.74360168269488,6.74545518923898,6.74566301077500,6.74421473994022,6.74110409841277,6.73632896499468,6.72989139035455,6.72179759838817,6.71205797420237,6.70068703877495,6.68770341039015,6.67312975299563,6.65699271167282,6.63932283545737,6.62015448779051,6.59952574492499,6.57747828265089,6.55405725174680,6.52931114260043,6.50329163947956,6.47605346496921,6.44765421512424,6.41815418591714,6.38761619158975,6.35610537554424,6.32368901443238,6.29043631612409,6.25641821225509,6.22170714606997,6.18637685629095,6.15050215775329,6.11415871955720,6.07742284149125,6.04037122948544,6.00308077085192,5.96562831006895,5.92809042585802,5.89054321029621,5.85306205069484,5.81572141496235,5.77859464115326,5.74175373188676,5.70526915429757,5.66920964615881,5.63364202879110,5.59863102734481,5.56423909901293,5.53052626970065,5.49754997964488,5.46536493844209,5.43402298990668,5.40357298714489,5.37406067819016,5.34552860250642,5.31801599862514,5.29155872314048,5.26618918124519,5.24193626894742,5.21882532706608,5.19687810705984,5.17611274870217,5.15654376957244,5.13818206629139,5.12103492738767,5.10510605764147,5.09039561371158,5.07690025081337,5.06461318017720,5.05352423698089,5.04361995841430,5.03488367150120,5.02729559027114,5.02083292184468,5.01546998096657,5.01117831249558,5.00792682133546,5.00568190926941,5.00440761814057,5.00406577880364,5.00461616525725,5.00601665335399,5.00822338347434,5.01119092654276,5.01487245275844,5.01921990241008,5.02418415814304,5.02971521804921,5.03576236895349,5.04227435927738,5.04919957086899,5.05648618919933,5.06408237133823,5.07193641113844,5.07999690107371,5.08821289019609,5.09653403769888,5.10491076159459,5.11329438204213,5.12163725888358,5.12989292297853,5.13801620095338,5.14596333301283,5.15369208349247,5.16116184386356,5.16833372793412,5.17517065902435,5.18163744892872,5.18770086851163,5.19332970981858,5.19849483961990,5.20316924433885,5.20732806635109,5.21094863167663,5.21401046911970,5.21649532094526,5.21838714521372,5.21967210992706,5.22033857917088,5.22037709146608,5.21978033057304,5.21854308901805,5.21666222463773,5.21413661046148,5.21096707827468,5.20715635622649,5.20270900086529,5.19763132400220,5.19193131481890,5.18561855764930,5.17870414587658,5.17120059239643,5.16312173710531,5.15448265187775,5.14529954350052,5.13558965503265,5.12537116606012,5.11466309231087,5.10348518509167,5.09185783100152,5.07980195236781,5.06733890884114,5.05449040057260,5.04127837338304,5.02772492631886,5.01385222197108,4.99968239991611,4.98523749361649,4.97053935109836,4.95560955970002,4.94046937516190,4.92513965530386,4.90964079851003,4.89399268721506,4.87821463655867,4.86232534834800,4.84634287043933,4.83028456162253,4.81416706206368,4.79800626933255,4.78181732001397,4.76561457687403,4.74941162152433,4.73322125250110,4.71705548864876,4.70092557767269,4.68484200970040,4.66881453566692,4.65285219031717,4.63696331959638,4.62115561217930,4.60543613486980,4.58981137158465,4.57428726561932,4.55886926487872,4.54356236974297,4.52837118322683,4.51329996308179,4.49835267548203,4.48353304992920,4.46884463500679,4.45429085461238,4.43987506429541,4.42560060732922,4.41147087014934,4.39748933679446,4.38365964199348,4.36998562255008,4.35647136668637,4.34312126101884,4.32994003485313,4.31693280149885,4.30410509632195,4.29146291126994,4.27901272562382,4.26676153275121,4.25471686265605,4.24288680014272,4.23127999843557,4.21990568811892,4.20877368128727,4.19789437082075,4.18727872472660,4.17693827551367,4.16688510459318,4.15713182172554,4.14769153955929,4.13857784333478,4.12980475585099,4.12138669781986,4.11333844375733,4.10567507358506,4.09841192014071,4.09156451281733,4.08514851757475,4.07917967358669,4.07367372680718,4.06864636075850,4.06411312485997,4.06008936063292,4.05659012613146,4.05363011896122,4.05122359825979,4.04938430602156,4.04812538815778,4.04745931568816,4.04739780646500,4.04795174783239,4.04913112062413,4.05094492490204,4.05340110783322,4.05650649409952,4.06026671922550,4.06468616620214,4.06976790577290,4.07551364073630,4.08192365460488,4.08899676494466,4.09673028170174,4.10511997080353,4.11416002330181,4.12384303030291,4.13415996390711,4.14510016435516,4.15665133355411,4.16879953512864,4.18152920111626,4.19482314539718,4.20866258392038,4.22302716175865,4.23789498699533,4.25324267141566,4.26904537794551,4.28527687475021,4.30190959587588,4.31891470828631,4.33626218511853,4.35392088495162,4.37185863685472,4.39004233095287,4.40843801422238,4.42701099120181,4.44572592928002,4.46454696819931,4.48343783338960,4.50236195272894,4.52128257630653,4.54016289874688,4.55896618363784,4.57765588959147,4.59619579745407,4.61455013817182,4.63268372081000,4.65056206021726,4.66815150382264,4.68541935705032,4.70233400683745,4.71886504274232,4.73498337513433,4.75066134996359,4.76587285961641,4.78059344937355,4.79480041900065,4.80847291901488,4.82159204118864,4.83414090286949,4.84610472471603,4.85747090147172,4.86822906542239,4.87837114220891,4.88789139869360,4.89678648260724,4.90505545373383,4.91269980642092,4.91972348323580];
r = r ./ r(end) .* 4.5;

%% Reverse
if reverse == true
    for i = 1:length(r)
        rr(i) = r(end-(i-1))
    end
    close all
    rr_throat = min(rr);
    throat = find(rr_throat == rr); %n-loaction of throat
    if length(throat)>1
        display("Warning: 2 locations for throat found")
    end
    throat = throat(1);
    zz_throat = z(throat)
    figure(1);
    plot([zz_throat zz_throat],[rr_throat -rr_throat],'--')
    hold on
    plot(z,rr,z,-rr)
    ylabel('Radius')
    xlabel('z')
    d_e_d_s = rr(end)/rr_throat
end