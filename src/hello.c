/**
@file hello.c
@brief Prints "Hello world !" in the uart terminal. Also does other examples 
*/

#include<uart.h>
#include "gpio.h"
#include "math.h"
#include <string.h>

/* Only enable one at a time, multiple enables will concede to the first defition */
#define BLINK_LED_EXAMPLE  	0
#define LED_COUNTER_EXAMPLE 1
#define DONUT_EXAMPLE		0

#if LED_COUNTER_EXAMPLE
	#define BIT_POS_0 	1 << 0
	#define BIT_POS_1 	1 << 1
	#define BIT_POS_2 	1 << 2
	#define BIT_POS_3 	1 << 3
#endif 

#if DONUT_EXAMPLE
/* Uncomment one of these to run individual examples */
// #define DONUT_FLOATING
// #define DONUT_LOOKUP
// #define DONUT_FIXED_PT
#define DONUT_SHORT_FIXED
#define M_PI		3.14159265358979323846
void print_donut(void);
#endif 

/** @fn void main()
 * @brief prints hello world and other demos
 */
int main()
{

	#if LED_COUNTER_EXAMPLE
		uint8_t counter = 1;
	#endif 

	#if BLINK_LED_EXAMPLE
		write_word(GPIO_DIRECTION_CNTRL_REG, GPIO16);
	#elif LED_COUNTER_EXAMPLE
		/* Set Blue LED's as outputs */
		write_word(GPIO_DIRECTION_CNTRL_REG, GPIO16 | GPIO19 |  GPIO22 |  GPIO25 );
		/* Clear any old data */
		write_word(GPIO_DATA_REG, 0X0); 
	#elif DONUT_EXAMPLE
		/* Change baud rate to 115200 for fast rendering */
		set_baud_rate(uart_instance[0], 115200);
	#endif

	while(1)
	{
		
			#if BLINK_LED_EXAMPLE
				/* Write to GPIO16 which is connected to LED*/
				write_word(GPIO_DATA_REG, GPIO16); 
				/* Terrible delay function, but it's roughly 1 second. Odd considering it's only looping a 1000 times in the func*/
				delay(1);
				/* Clear all*/
				write_word(GPIO_DATA_REG, 0X0); 
				delay(1);
			#elif LED_COUNTER_EXAMPLE

				if(counter & BIT_POS_0)
				{
					write_word(GPIO_DATA_REG, read_word(GPIO_DATA_REG) | GPIO16); 
				}
				else
				{
					write_word(GPIO_DATA_REG, (read_word(GPIO_DATA_REG) & ~GPIO16)); 
				}

				if(counter & BIT_POS_1)
				{
					write_word(GPIO_DATA_REG, read_word(GPIO_DATA_REG) | GPIO19); 
				}
				else
				{
					write_word(GPIO_DATA_REG, (read_word(GPIO_DATA_REG) & ~GPIO19)); 
				}

				if(counter & BIT_POS_2)
				{
					write_word(GPIO_DATA_REG, read_word(GPIO_DATA_REG) | GPIO22); 
				}
				else
				{
					write_word(GPIO_DATA_REG, read_word(GPIO_DATA_REG) & ~GPIO22); 
				}

				if(counter & BIT_POS_3)
				{
					write_word(GPIO_DATA_REG, read_word(GPIO_DATA_REG) | GPIO25); 
				}
				else
				{
					write_word(GPIO_DATA_REG, read_word(GPIO_DATA_REG) & ~GPIO25); 
				}

				counter++;
				if(counter > 15)
				{	
					/* Wait just a little longer when count is done */
					delay_loop(1000,1000);
					counter = 0;
				}
				else
				{
					/* Roughly 0.4 seconds */
					delay_loop(1000,400);
					// printf("Test");
					// printf("\x1b[H");

				}
			#elif DONUT_EXAMPLE
				print_donut();
			#endif

	}
}


#if DONUT_EXAMPLE

/* Readable version Floating point  donut */
void donut_readable_float() {
	int screen_width = 40, screen_height=40;
	uint32_t i;
	float A =0, B=0;
	const float theta_spacing = 0.3;
	const float phi_spacing   = 0.1;

	const float R1 = 1;
	const float R2 = 2;
	const float K2 = 5;
	// Calculate K1 based on screen size: the maximum x-distance occurs
	// roughly at the edge of the torus, which is at x=R1+R2, z=0.  we
	// want that to be displaced 3/8ths of the width of the screen, which
	// is 3/4th of the way from the center to the side of the screen.
	// screen_width*3/8 = K1*(R1+R2)/(K2+0)
	// screen_width*K2*3/(8*(R1+R2)) = K1
	const float K1 = screen_width*K2*3/(8*(R1+R2));
	char output[screen_width][screen_height];
	float zbuffer[screen_width][screen_height];


	while(1)
	{
		memset(output,32,screen_width*screen_height);
		memset(zbuffer,0,screen_width*screen_height*4);

		// precompute sines and cosines of A and B
		float cosA = cos(A), sinA = sin(A);
		float cosB = cos(B), sinB = sin(B);


		// theta goes around the cross-sectional circle of a torus
		for (float theta=0; theta < 2*M_PI; theta += theta_spacing) {
			// precompute sines and cosines of theta
			float costheta = cos(theta), sintheta = sin(theta);

			// phi goes around the center of revolution of a torus
			for(float phi=0; phi < 2*M_PI; phi += phi_spacing) {
			// precompute sines and cosines of phi
			float cosphi = cos(phi), sinphi = sin(phi);
			
			// the x,y coordinate of the circle, before revolving (factored
			// out of the above equations)
			float circlex = R2 + R1*costheta;
			float circley = R1*sintheta;

			// final 3D (x,y,z) coordinate after rotations, directly from
			// our math above
			float x = circlex*(cosB*cosphi + sinA*sinB*sinphi)
				- circley*cosA*sinB; 
			float y = circlex*(sinB*cosphi - sinA*cosB*sinphi)
				+ circley*cosA*cosB;
			float z = K2 + cosA*circlex*sinphi + circley*sinA;
			float ooz = 1/z;  // "one over z"
			
			// x and y projection.  note that y is negated here, because y
			// goes up in 3D space but down on 2D displays.
			int xp = (int) (screen_width/2 + K1*ooz*x);
			int yp = (int) (screen_height/2 - K1*ooz*y);
			
			// calculate luminance.  ugly, but correct.
			float L = cosphi*costheta*sinB - cosA*costheta*sinphi -
				sinA*sintheta + cosB*(cosA*sintheta - costheta*sinA*sinphi);
			// L ranges from -sqrt(2) to +sqrt(2).  If it's < 0, the surface
			// is pointing away from us, so we won't bother trying to plot it.
			if (L > 0) {
				// test against the z-buffer.  larger 1/z means the pixel is
				// closer to the viewer than what's already plotted.
				if(ooz > zbuffer[xp][yp]) {
				zbuffer[xp][yp] = ooz;
				int luminance_index = L*8;
				// luminance_index is now in the range 0..11 (8*sqrt(2) = 11.3)
				// now we lookup the character corresponding to the
				// luminance and plot it in our output:
				output[xp][yp] = ".,-~:;=!*#$@"[luminance_index];
				}
			}
			}
		}

		// now, dump output[] to the screen.
		// bring cursor to "home" location, in just about any currently-used
		// terminal emulation mode
		printf("\x1b[H");
		for (int j = 0; j < screen_height; j++) {

			for (int i = 0; i < screen_width; i++) {
			putchar(output[i][j]);
			}
			putchar('\n');
		}
			A+=0.10;
			B+=0.10;
	}
	
}



/* Lookup table for sine and cosine */
float sin_lookup[315] = {0,0.019998667,0.039989334,0.059964006,0.079914694,0.099833417,0.119712207,0.139543115,0.159318207,0.179029573,0.198669331,0.218229623,0.237702626,0.257080552,0.276355649,0.295520207,0.314566561,0.333487092,0.352274233,0.370920469,0.389418342,0.407760453,0.425939465,0.443948107,0.461779176,0.479425539,0.496880138,0.514135992,0.531186198,0.548023937,0.564642473,0.581035161,0.597195441,0.613116852,0.628793024,0.644217687,0.659384672,0.674287912,0.688921445,0.703279419,0.717356091,0.73114583,0.74464312,0.757842563,0.770738879,0.78332691,0.79560162,0.8075581,0.819191568,0.83049737,0.841470985,0.852108022,0.862404227,0.872355482,0.881957807,0.89120736,0.900100442,0.908633496,0.916803109,0.924606012,0.932039086,0.939099356,0.945783999,0.952090342,0.95801586,0.963558185,0.9687151,0.973484542,0.977864602,0.98185353,0.98544973,0.988651763,0.991458348,0.993868363,0.995880845,0.997494987,0.998710144,0.999525831,0.99994172,0.999957646,0.999573603,0.998789743,0.997606381,0.99602399,0.994043202,0.99166481,0.988889766,0.985719179,0.982154317,0.978196607,0.973847631,0.969109129,0.963982996,0.958471283,0.952576194,0.946300088,0.939645474,0.932615014,0.925211521,0.917437955,0.909297427,0.900793192,0.891928651,0.882707351,0.87313298,0.863209367,0.852940482,0.842330432,0.831383461,0.820103948,0.808496404,0.796565472,0.784315925,0.771752662,0.758880708,0.745705212,0.732231444,0.718464793,0.704410766,0.690074984,0.675463181,0.660581201,0.645434998,0.63003063,0.614374258,0.598472144,0.58233065,0.56595623,0.549355436,0.532534908,0.515501372,0.498261642,0.480822615,0.463191265,0.445374645,0.42737988,0.40921417,0.390884779,0.372399039,0.353764345,0.33498815,0.316077964,0.297041351,0.277885926,0.25861935,0.239249329,0.219783612,0.200229985,0.180596268,0.160890315,0.141120008,0.121293255,0.101417986,0.081502152,0.061553717,0.041580662,0.021590976,0.001592653,-0.018406307,-0.038397905,-0.058374143,-0.078327033,-0.098248594,-0.118130856,-0.137965867,-0.157745694,-0.177462425,-0.197108173,-0.21667508,-0.236155321,-0.255541102,-0.27482467,-0.293998312,-0.313054359,-0.331985188,-0.350783228,-0.369440959,-0.387950918,-0.406305702,-0.424497969,-0.442520443,-0.460365915,-0.478027246,-0.495497373,-0.512769307,-0.529836141,-0.546691047,-0.563327284,-0.579738198,-0.595917224,-0.611857891,-0.627553823,-0.642998742,-0.65818647,-0.673110932,-0.687766159,-0.702146289,-0.716245569,-0.730058361,-0.743579139,-0.756802495,-0.769723141,-0.782335907,-0.79463575,-0.806617749,-0.818277111,-0.829609174,-0.840609404,-0.851273401,-0.8615969,-0.871575772,-0.881206026,-0.890483809,-0.89940541,-0.907967261,-0.916165937,-0.923998159,-0.931460794,-0.938550857,-0.945265512,-0.951602074,-0.957558007,-0.963130931,-0.968318614,-0.973118983,-0.977530118,-0.981550253,-0.985177782,-0.988411252,-0.991249371,-0.993691004,-0.995735173,-0.997381062,-0.998628011,-0.999475523,-0.999923258,-0.999971036,-0.99961884,-0.998866809,-0.997715246,-0.996164609,-0.99421552,-0.991868757,-0.989125261,-0.985986127,-0.982452613,-0.97852613,-0.97420825,-0.969500699,-0.964405362,-0.958924275,-0.953059631,-0.946813776,-0.940189208,-0.933188576,-0.925814682,-0.918070475,-0.909959051,-0.901483656,-0.892647679,-0.883454656,-0.873908262,-0.864012316,-0.853770778,-0.843187742,-0.832267442,-0.821014247,-0.809432656,-0.797527304,-0.785302951,-0.772764488,-0.759916929,-0.746765413,-0.733315201,-0.719571673,-0.705540326,-0.691226772,-0.676636736,-0.661776055,-0.646650672,-0.631266638,-0.615630105,-0.599747329,-0.583624661,-0.567268552,-0.550685543,-0.533882266,-0.516865444,-0.499641883,-0.482218472,-0.464602179,-0.446800052,-0.428819211,-0.410666848,-0.392350224,-0.373876665,-0.35525356,-0.336488358,-0.317588566,-0.298561742,-0.279415498,-0.260157491,-0.240795425,-0.221337044,-0.201790131,-0.182162504,-0.162462015,-0.142696544,-0.122873995,-0.103002299,-0.083089403,-0.063143272,-0.043171885,-0.02318323,-0.003185302};
float cos_lookup[315] = {1,0.999800007,0.999200107,0.99820054,0.996801706,0.995004165,0.992808636,0.990215996,0.987227283,0.983843693,0.980066578,0.975897449,0.971337975,0.966389978,0.961055438,0.955336489,0.949235418,0.942754666,0.935896824,0.928664636,0.921060994,0.91308894,0.904751663,0.896052498,0.886994923,0.877582562,0.86781918,0.857708681,0.847255111,0.83646265,0.825335615,0.813878457,0.802095758,0.789992231,0.777572719,0.764842187,0.751805729,0.738468559,0.724836011,0.710913538,0.696706709,0.682221207,0.667462826,0.652437468,0.637151144,0.621609968,0.605820157,0.589788025,0.573519986,0.557022547,0.540302306,0.523365951,0.506220257,0.488872082,0.471328364,0.453596121,0.435682446,0.417594504,0.399339529,0.380924824,0.362357754,0.343645746,0.324796284,0.305816908,0.28671521,0.267498829,0.248175452,0.228752808,0.209238666,0.189640831,0.169967143,0.15022547,0.130423709,0.11056978,0.090671624,0.070737202,0.050774485,0.030791459,0.010796117,-0.009203543,-0.029199522,-0.049183822,-0.069148449,-0.089085417,-0.108986752,-0.128844494,-0.1486507,-0.168397448,-0.188076839,-0.207681002,-0.227202095,-0.24663231,-0.265963876,-0.285189059,-0.304300171,-0.323289567,-0.342149651,-0.36087288,-0.379451765,-0.397878874,-0.416146837,-0.434248346,-0.452176162,-0.469923114,-0.487482102,-0.504846105,-0.522008175,-0.538961449,-0.555699146,-0.572214571,-0.588501117,-0.604552271,-0.620361612,-0.635922817,-0.651229661,-0.666276021,-0.681055881,-0.695563326,-0.709792556,-0.723737879,-0.737393716,-0.750754605,-0.763815202,-0.776570284,-0.789014747,-0.801143616,-0.812952037,-0.824435289,-0.835588777,-0.846408041,-0.856888753,-0.867026721,-0.87681789,-0.886258344,-0.895344306,-0.904072142,-0.912438361,-0.920439618,-0.92807271,-0.935334586,-0.942222341,-0.948733219,-0.954864616,-0.960614081,-0.965979312,-0.970958165,-0.975548648,-0.979748924,-0.983557313,-0.986972293,-0.989992497,-0.992616717,-0.994843903,-0.996673166,-0.998103772,-0.99913515,-0.999766888,-0.999998732,-0.99983059,-0.999262529,-0.998294776,-0.996927718,-0.995161903,-0.992998037,-0.990436984,-0.98747977,-0.984127577,-0.980381746,-0.976243776,-0.971715321,-0.966798193,-0.961494358,-0.955805939,-0.94973521,-0.943284599,-0.936456687,-0.929254205,-0.921680034,-0.913737203,-0.905428889,-0.896758416,-0.887729252,-0.878345007,-0.868609437,-0.858526434,-0.848100032,-0.837334401,-0.826233848,-0.814802812,-0.803045866,-0.790967712,-0.778573182,-0.765867232,-0.752854947,-0.739541529,-0.725932304,-0.712032716,-0.697848325,-0.683384804,-0.668647937,-0.653643621,-0.638377856,-0.622856748,-0.607086506,-0.591073437,-0.574823947,-0.558344534,-0.541641792,-0.5247224,-0.507593126,-0.490260821,-0.472732419,-0.45501493,-0.437115441,-0.419041112,-0.400799172,-0.382396918,-0.36384171,-0.34514097,-0.326302178,-0.30733287,-0.288240633,-0.269033103,-0.249717964,-0.230302941,-0.210795799,-0.191204343,-0.171536407,-0.151799858,-0.132002592,-0.112152527,-0.092257602,-0.072325775,-0.052365019,-0.032383318,-0.012388663,0.007610946,0.027607511,0.047593034,0.06755952,0.087498983,0.107403448,0.127264953,0.147075554,0.166827326,0.186512369,0.206122811,0.225650805,0.245088543,0.264428248,0.283662185,0.302782662,0.321782029,0.340652688,0.35938709,0.377977743,0.396417209,0.414698114,0.432813144,0.450755056,0.468516671,0.486090886,0.503470671,0.520649075,0.537619226,0.554374336,0.570907704,0.587212717,0.603282852,0.619111682,0.634692876,0.650020201,0.665087527,0.679888826,0.694418179,0.708669774,0.722637911,0.736317002,0.749701576,0.762786279,0.775565879,0.788035262,0.800189441,0.812023555,0.823532871,0.834712785,0.845558824,0.856066652,0.866232064,0.876050995,0.885519517,0.894633843,0.903390328,0.911785468,0.919815906,0.927478431,0.934769976,0.941687626,0.948228613,0.954390322,0.960170287,0.965566196,0.970575893,0.975197371,0.979428784,0.983268438,0.986714799,0.989766486,0.99242228,0.994681118,0.996542097,0.998004473,0.99906766,0.999731233,0.999994927};

#define SCREEN_WIDTH 40u
#define SCREEN_HEIGHT 40u

/* Donut with Sine lookup point math */
void donut_readable_fixed_lookup() {
	
	int32_t A =0, B=0;

	const int32_t theta_spacing = 15;
	const int32_t phi_spacing   = 5;

	const float R1 = 1;
	const float R2 = 2;
	const float K2 = 5;
	// Calculate K1 based on screen size: the maximum x-distance occurs
	// roughly at the edge of the torus, which is at x=R1+R2, z=0.  we
	// want that to be displaced 3/8ths of the width of the screen, which
	// is 3/4th of the way from the center to the side of the screen.
	// screen_width*3/8 = K1*(R1+R2)/(K2+0)
	// screen_width*K2*3/(8*(R1+R2)) = K1
	const float K1 = SCREEN_WIDTH*K2*3/(8*(R1+R2));
	
	/* Precalculating some R1,R2 and K2 for readability */ 
	// const int K1 = q_mul(SCREEN_WIDTH, (int32_t) 0.625*MUL_QTY);

	char output[SCREEN_WIDTH][SCREEN_HEIGHT];
	float zbuffer[SCREEN_WIDTH][SCREEN_HEIGHT];

	while(1)
	{
		memset(output,32,SCREEN_WIDTH*SCREEN_HEIGHT);
		memset(zbuffer,0,SCREEN_WIDTH*SCREEN_HEIGHT*4);

		// precompute sines and cosines of A and B
		float cosA = cos_lookup[A], sinA = sin_lookup[A];
		float cosB = cos_lookup[B], sinB = sin_lookup[B];


		// theta goes around the cross-sectional circle of a torus
		for (uint32_t theta=0; theta < 315; theta += theta_spacing) {
			// precompute sines and cosines of theta
			float costheta = cos_lookup[theta], sintheta = sin_lookup[theta];

			// phi goes around the center of revolution of a torus
			for(uint32_t phi=0; phi < 315; phi += phi_spacing) {
			// precompute sines and cosines of phi
			float cosphi = cos_lookup[phi], sinphi = sin_lookup[phi];
			
			// the x,y coordinate of the circle, before revolving (factored
			// out of the above equations)
			float circlex = R2 + R1*costheta;
			float circley = R1*sintheta;

			// final 3D (x,y,z) coordinate after rotations, directly from
			// our math above
			float x = circlex*(cosB*cosphi + sinA*sinB*sinphi)
				- circley*cosA*sinB; 
			float y = circlex*(sinB*cosphi - sinA*cosB*sinphi)
				+ circley*cosA*cosB;
			float z = K2 + cosA*circlex*sinphi + circley*sinA;
			float ooz = 1/z;  // "one over z"

			// x and y projection.  note that y is negated here, because y
			// goes up in 3D space but down on 2D displays.
			int xp = (int) (SCREEN_WIDTH/2 + K1*ooz*x);
			int yp = (int) (SCREEN_HEIGHT/2 - K1*ooz*y);
			
			// calculate luminance.  ugly, but correct.
			float L = cosphi*costheta*sinB - cosA*costheta*sinphi -
				sinA*sintheta + cosB*(cosA*sintheta - costheta*sinA*sinphi);
			// L ranges from -sqrt(2) to +sqrt(2).  If it's < 0, the surface
			// is pointing away from us, so we won't bother trying to plot it.
			if (L > 0) {
				// test against the z-buffer.  larger 1/z means the pixel is
				// closer to the viewer than what's already plotted.
				if(ooz > zbuffer[xp][yp]) {
				zbuffer[xp][yp] = ooz;
				int luminance_index = L*8;
				// luminance_index is now in the range 0..11 (8*sqrt(2) = 11.3)
				// now we lookup the character corresponding to the
				// luminance and plot it in our output:
				output[xp][yp] = ".,-~:;=!*#$@"[luminance_index];
				}
			}

			}
		}

		// now, dump output[] to the screen.
		// bring cursor to "home" location, in just about any currently-used
		// terminal emulation mode
		printf("\x1b[H");
		for (int j = 0; j < SCREEN_HEIGHT; j++) {
			for (int i = 0; i < SCREEN_WIDTH; i++) {
			putchar(output[i][j]);
			}
			putchar('\n');
		}
			A+=5;
			/* Check for rollovers to take care of array lookups */
			if(A > 314)
			{
				A -= 314;
			}
			B+=5;
			/* Check for rollovers to take care of array lookups */
			if(B > 314)
			{
				B -= 314;
			}
	}
	
}


/* Q15.16 */
#define SHIFT_QTY 16
#define MUL_QTY (1 << SHIFT_QTY)

int32_t sin_lookup_fixed[315];
int32_t cos_lookup_fixed[315];

int32_t q_mul(int32_t t1, int32_t t2)
{
	int64_t acc = 0;
	acc = (int64_t) t1*(int64_t) t2;
	return((int32_t) (acc >> SHIFT_QTY));
}

int32_t q_div(int32_t t1, int32_t t2)
{
	int64_t acc = 0;
	
	acc =  (int64_t) t1*MUL_QTY;

	return((int32_t) (acc/t2));
}

/* Helpful if you want to see output values per operation */
// #define DEBUG_REFERENCE_ENABLE

/* Donut with 32-bit fixed point math */
void donut_readable_all_fixed() {
	int32_t A =0, B=0;
	const int32_t theta_spacing = 15;
	const int32_t phi_spacing   = 5;
	uint32_t i = 0;
	const int32_t R1 = 1<<SHIFT_QTY;
	const int32_t R2 = 2<<SHIFT_QTY;
	const int32_t K2 = 5<<SHIFT_QTY;
	#ifdef DEBUG_REFERENCE_ENABLE
	const float R1_ref = 1;
	const float R2_ref = 2;
	const float K2_ref = 5;
	#endif
	// Calculate K1 based on screen size: the maximum x-distance occurs
	// roughly at the edge of the torus, which is at x=R1+R2, z=0.  we
	// want that to be displaced 3/8ths of the width of the screen, which
	// is 3/4th of the way from the center to the side of the screen.
	// screen_width*3/8 = K1*(R1+R2)/(K2+0)
	// screen_width*K2*3/(8*(R1+R2)) = K1
	// const float K1 = screen_width*K2*3/(8*(R1+R2));
	#ifdef DEBUG_REFERENCE_ENABLE
	const float K1_ref = SCREEN_WIDTH*K2_ref*3/(8*(R1_ref+R2_ref));
	#endif

	const uint32_t K1 = q_mul(SCREEN_WIDTH<<SHIFT_QTY, (0.625f*MUL_QTY));

	char output[SCREEN_WIDTH][SCREEN_HEIGHT];
	int32_t zbuffer[SCREEN_WIDTH][SCREEN_HEIGHT];

	/* Precalc lookups */
	for(i=0;i < 315; i++)
	{
		sin_lookup_fixed[i] = (int32_t) (MUL_QTY*sin_lookup[i]);
		cos_lookup_fixed[i] = (int32_t) (MUL_QTY*cos_lookup[i]);

	}

	while(1)
	{
		memset(output,32,SCREEN_WIDTH*SCREEN_HEIGHT);
		memset(zbuffer,0,SCREEN_WIDTH*SCREEN_HEIGHT*4);

		// precompute sines and cosines of A and B
		int32_t cosA = cos_lookup_fixed[A], sinA = sin_lookup_fixed[A];
		int32_t cosB = cos_lookup_fixed[B], sinB = sin_lookup_fixed[B];
		#ifdef DEBUG_REFERENCE_ENABLE
		float cosA_ref = cos_lookup[A], sinA_ref = sin_lookup[A];
		float cosB_ref = cos_lookup[B], sinB_ref = sin_lookup[B];
		// printf("%f, %f,%ld, %ld \n", cosA_ref,sinA_ref,cosA,sinA);
		// printf("%f, %f,%ld, %ld \n", cosB_ref,sinB_ref,cosB,sinB);
		#endif


		
		// theta goes around the cross-sectional circle of a torus
		for (uint32_t theta=0; theta < 315; theta += theta_spacing) {
			// precompute sines and cosines of theta
			int32_t costheta = cos_lookup_fixed[theta], sintheta = sin_lookup_fixed[theta];
			// phi goes around the center of revolution of a torus
			for(uint32_t phi=0; phi < 315; phi += phi_spacing) {
			// precompute sines and cosines of phi
			int32_t cosphi = cos_lookup_fixed[phi], sinphi = sin_lookup_fixed[phi];
			#ifdef DEBUG_REFERENCE_ENABLE
			float costheta_ref = cos_lookup[theta], sintheta_ref = sin_lookup[theta];
			float cosphi_ref = cos_lookup[phi], sinphi_ref = sin_lookup[phi];
			// printf("%f, %f, %f, %f, %ld, %ld, %ld, %ld \n", costheta_ref,sintheta_ref,cosphi_ref,sinphi_ref,costheta,sintheta,cosphi,sinphi);
			#endif


			// the x,y coordinate of the circle, before revolving (factored
			// out of the above equations)
			int32_t circlex = R2 + q_mul(R1,costheta);
			int32_t circley = q_mul(R1,sintheta);
			#ifdef DEBUG_REFERENCE_ENABLE
			float circlex_ref = R2_ref + R1_ref*costheta_ref;
			float circley_ref = R1_ref*sintheta_ref;
			#endif
			// printf("%f, %f,%ld, %ld \n", circlex_ref,circley_ref,circlex,circley);

			// final 3D (x,y,z) coordinate after rotations, directly from
			// our math above
			// float x = circlex*(cosB*cosphi + sinA*sinB*sinphi)
			// 	- circley*cosA*sinB; 

			int32_t x = q_mul(circlex,(q_mul(cosB,cosphi) + q_mul(q_mul(sinA,sinB),sinphi)))
				- q_mul(circley,q_mul(cosA,sinB)); 
			int32_t y = q_mul(circlex,(q_mul(sinB,cosphi) - q_mul(q_mul(sinA,cosB),sinphi)))
				+ q_mul(circley,q_mul(cosA,cosB));
			int32_t z = K2 + q_mul(q_mul(cosA,circlex),sinphi) + q_mul(circley,sinA);
			int32_t ooz = q_div((1<<SHIFT_QTY),z);  // "one over z"

			#ifdef DEBUG_REFERENCE_ENABLE
			float x_ref = circlex_ref*(cosB_ref*cosphi_ref + sinA_ref*sinB_ref*sinphi_ref)
				- circley_ref*cosA_ref*sinB_ref; 
			float y_ref = circlex_ref*(sinB_ref*cosphi_ref - sinA_ref*cosB_ref*sinphi_ref)
				+ circley_ref*cosA_ref*cosB_ref;
			float z_ref = K2_ref + cosA_ref*circlex_ref*sinphi_ref + circley_ref*sinA_ref;
			float ooz_ref = 1/z_ref;  // "one over z"
			// printf("%f, %f,%ld, %ld \n", z_ref,ooz_ref,z,ooz);
			#endif
			// x and y projection.  note that y is negated here, because y
			// goes up in 3D space but down on 2D displays.
			int32_t xp = (int) ((SCREEN_WIDTH<<(SHIFT_QTY-1)) + q_mul(q_mul(K1,ooz),x))>>SHIFT_QTY;
			int32_t yp = (int) ((SCREEN_HEIGHT<<(SHIFT_QTY-1)) - q_mul(q_mul(K1,ooz),y))>>SHIFT_QTY;
			
			#ifdef DEBUG_REFERENCE_ENABLE
			int xp_ref = (int) (SCREEN_WIDTH/2 + K1_ref*ooz_ref*x_ref);
			int yp_ref = (int) (SCREEN_HEIGHT/2 - K1_ref*ooz_ref*y_ref);
			// printf("%ld, %ld,%ld, %ld \n", xp_ref,yp_ref,xp,yp);
			#endif


			// calculate luminance.  ugly, but correct.
			// float L = cosphi*costheta*sinB - cosA*costheta*sinphi -
			// 	sinA*sintheta + cosB*(cosA*sintheta - costheta*sinA*sinphi);

			int32_t L = q_mul(q_mul(cosphi,costheta),sinB) - q_mul(q_mul(cosA,costheta),sinphi) -
				q_mul(sinA,sintheta) + q_mul(cosB,q_mul(cosA,sintheta) - q_mul(q_mul(costheta,sinA),sinphi));

			#ifdef DEBUG_REFERENCE_ENABLE
			float L_ref = cosphi_ref*costheta_ref*sinB_ref - cosA_ref*costheta_ref*sinphi_ref -
				sinA_ref*sintheta_ref + cosB_ref*(cosA_ref*sintheta_ref - costheta_ref*sinA_ref*sinphi_ref);
			// printf("%f, %ld \n", L_ref,L);
			#endif
			// L ranges from -sqrt(2) to +sqrt(2).  If it's < 0, the surface
			// is pointing away from us, so we won't bother trying to plot it.
			if (L > 0) {
				// test against the z-buffer.  larger 1/z means the pixel is
				// closer to the viewer than what's already plotted.
				if(ooz > zbuffer[xp][yp]) {
				zbuffer[xp][yp] = ooz;
				/* Here, the multiply by 8 has been transformed by removing some shift  */
				L >>= (SHIFT_QTY-3);
				int luminance_index = L;

				// luminance_index is now in the range 0..11 (8*sqrt(2) = 11.3)
				// now we lookup the character corresponding to the
				// luminance and plot it in our output:
				output[xp][yp] = ".,-~:;=!*#$@"[luminance_index];
				}
			}
			}
		}

		// now, dump output[] to the screen.
		// bring cursor to "home" location, in just about any currently-used
		// terminal emulation mode
		printf("\x1b[H");
		for (int j = 0; j < SCREEN_HEIGHT; j++) {
			for (int i = 0; i < SCREEN_WIDTH; i++) {
			putchar(output[i][j]);
			}
			putchar('\n');
		}
		A+=5;
		/* Check for rollovers to take care of array lookups */
		if(A > 314)
		{
			A -= 314;
		}
		B+=5;
		/* Check for rollovers to take care of array lookups */
		if(B > 314)
		{
			B -= 314;
		}
	}
	
}




/* Q7.8 */
#define SHIFT_QTY_U16 8
#define MUL_QTY_U16 (1 << SHIFT_QTY_U16)

int16_t sin_lookup_fixed_u16[315];
int16_t cos_lookup_fixed_u16[315];

int16_t q_mul_u16(int16_t t1, int16_t t2)
{
	int32_t acc = 0;
	acc =t1*t2;
	return((int16_t) (acc >> SHIFT_QTY_U16));
}

int16_t q_div_u16(int16_t t1, int16_t t2)
{
	int32_t acc = 0;
	
	acc =  (int32_t) t1*MUL_QTY_U16;

	return((int16_t) (acc/t2));
}

/* Donut using int16 math only */
/* Also improves the printing algorithm by reversing height/width for easier access and a direct printf */
void donut_readable_all_fixed_u16() {
	int16_t A =0, B=0;
	const int16_t theta_spacing = 15;
	const int16_t phi_spacing   = 5;
	int16_t i = 0;
	const int16_t R1 = 1<<SHIFT_QTY_U16;
	const int16_t R2 = 2<<SHIFT_QTY_U16;
	const int16_t K2 = 5<<SHIFT_QTY_U16;
	#ifdef DEBUG_REFERENCE_ENABLE
	const float R1_ref = 1;
	const float R2_ref = 2;
	const float K2_ref = 5;
	#endif
	// Calculate K1 based on screen size: the maximum x-distance occurs
	// roughly at the edge of the torus, which is at x=R1+R2, z=0.  we
	// want that to be displaced 3/8ths of the width of the screen, which
	// is 3/4th of the way from the center to the side of the screen.
	// screen_width*3/8 = K1*(R1+R2)/(K2+0)
	// screen_width*K2*3/(8*(R1+R2)) = K1
	// const float K1 = screen_width*K2*3/(8*(R1+R2));
	#ifdef DEBUG_REFERENCE_ENABLE
	const float K1_ref = SCREEN_WIDTH*K2_ref*3/(8*(R1_ref+R2_ref));
	#endif

	const int16_t K1 = q_mul_u16(SCREEN_WIDTH<<SHIFT_QTY_U16, (0.625f*MUL_QTY_U16));

	char output[SCREEN_HEIGHT][SCREEN_WIDTH];
	int16_t zbuffer[SCREEN_HEIGHT][SCREEN_WIDTH];

	/* Precalc lookups */
	for(i=0;i < 315; i++)
	{
		sin_lookup_fixed_u16[i] = (int16_t) (MUL_QTY_U16*sin_lookup[i]);
		cos_lookup_fixed_u16[i] = (int16_t) (MUL_QTY_U16*cos_lookup[i]);
	}

	while(1)
	{
		memset(output,32,SCREEN_WIDTH*SCREEN_HEIGHT);
		memset(zbuffer,0,SCREEN_WIDTH*SCREEN_HEIGHT*2);

		// precompute sines and cosines of A and B
		int16_t cosA = cos_lookup_fixed_u16[A], sinA = sin_lookup_fixed_u16[A];
		int16_t cosB = cos_lookup_fixed_u16[B], sinB = sin_lookup_fixed_u16[B];
		#ifdef DEBUG_REFERENCE_ENABLE
		float cosA_ref = cos_lookup[A], sinA_ref = sin_lookup[A];
		float cosB_ref = cos_lookup[B], sinB_ref = sin_lookup[B];
		// printf("%f, %f,%ld, %ld \n", cosA_ref,sinA_ref,cosA,sinA);
		// printf("%f, %f,%ld, %ld \n", cosB_ref,sinB_ref,cosB,sinB);
		#endif

		
		// theta goes around the cross-sectional circle of a torus
		for (int16_t theta=0; theta < 315; theta += theta_spacing) {
			// precompute sines and cosines of theta
			int16_t costheta = cos_lookup_fixed_u16[theta], sintheta = sin_lookup_fixed_u16[theta];
			// phi goes around the center of revolution of a torus
			for(int16_t phi=0; phi < 315; phi += phi_spacing) {
			// precompute sines and cosines of phi
			int16_t cosphi = cos_lookup_fixed_u16[phi], sinphi = sin_lookup_fixed_u16[phi];
			#ifdef DEBUG_REFERENCE_ENABLE
			float costheta_ref = cos_lookup[theta], sintheta_ref = sin_lookup[theta];
			float cosphi_ref = cos_lookup[phi], sinphi_ref = sin_lookup[phi];
			// printf("%f, %f, %f, %f, %ld, %ld, %ld, %ld \n", costheta_ref,sintheta_ref,cosphi_ref,sinphi_ref,costheta,sintheta,cosphi,sinphi);
			#endif


			// the x,y coordinate of the circle, before revolving (factored
			// out of the above equations)
			int16_t circlex = R2 + q_mul_u16(R1,costheta);
			int16_t circley = q_mul_u16(R1,sintheta);
			#ifdef DEBUG_REFERENCE_ENABLE
			float circlex_ref = R2_ref + R1_ref*costheta_ref;
			float circley_ref = R1_ref*sintheta_ref;
			#endif
			// printf("%f, %f,%ld, %ld \n", circlex_ref,circley_ref,circlex,circley);

			// final 3D (x,y,z) coordinate after rotations, directly from
			// our math above
			// float x = circlex*(cosB*cosphi + sinA*sinB*sinphi)
			// 	- circley*cosA*sinB; 

			int16_t x = q_mul_u16(circlex,(q_mul_u16(cosB,cosphi) + q_mul_u16(q_mul_u16(sinA,sinB),sinphi)))
				- q_mul_u16(circley,q_mul_u16(cosA,sinB)); 
			int16_t y = q_mul_u16(circlex,(q_mul_u16(sinB,cosphi) - q_mul_u16(q_mul_u16(sinA,cosB),sinphi)))
				+ q_mul_u16(circley,q_mul_u16(cosA,cosB));
			int16_t z = K2 + q_mul_u16(q_mul_u16(cosA,circlex),sinphi) + q_mul_u16(circley,sinA);
			int16_t ooz = q_div_u16((1<<SHIFT_QTY_U16),z);  // "one over z"

			#ifdef DEBUG_REFERENCE_ENABLE
			float x_ref = circlex_ref*(cosB_ref*cosphi_ref + sinA_ref*sinB_ref*sinphi_ref)
				- circley_ref*cosA_ref*sinB_ref; 
			float y_ref = circlex_ref*(sinB_ref*cosphi_ref - sinA_ref*cosB_ref*sinphi_ref)
				+ circley_ref*cosA_ref*cosB_ref;
			float z_ref = K2_ref + cosA_ref*circlex_ref*sinphi_ref + circley_ref*sinA_ref;
			float ooz_ref = 1/z_ref;  // "one over z"
			// printf("%f, %f,%ld, %ld \n", z_ref,ooz_ref,z,ooz);
			#endif
			// x and y projection.  note that y is negated here, because y
			// goes up in 3D space but down on 2D displays.
			int16_t xp = (int16_t) ((SCREEN_WIDTH<<(SHIFT_QTY_U16-1)) + q_mul_u16(q_mul_u16(K1,ooz),x))>>SHIFT_QTY_U16;
			int16_t yp = (int16_t) ((SCREEN_HEIGHT<<(SHIFT_QTY_U16-1)) - q_mul_u16(q_mul_u16(K1,ooz),y))>>SHIFT_QTY_U16;
			
			#ifdef DEBUG_REFERENCE_ENABLE
			int xp_ref = (int) (SCREEN_WIDTH/2 + K1_ref*ooz_ref*x_ref);
			int yp_ref = (int) (SCREEN_HEIGHT/2 - K1_ref*ooz_ref*y_ref);
			// printf("%ld, %ld,%ld, %ld \n", xp_ref,yp_ref,xp,yp);
			#endif


			// calculate luminance.  ugly, but correct.
			// float L = cosphi*costheta*sinB - cosA*costheta*sinphi -
			// 	sinA*sintheta + cosB*(cosA*sintheta - costheta*sinA*sinphi);

			int16_t L = q_mul_u16(q_mul_u16(cosphi,costheta),sinB) - q_mul_u16(q_mul_u16(cosA,costheta),sinphi) -
				q_mul_u16(sinA,sintheta) + q_mul_u16(cosB,q_mul_u16(cosA,sintheta) - q_mul_u16(q_mul_u16(costheta,sinA),sinphi));

			#ifdef DEBUG_REFERENCE_ENABLE
			float L_ref = cosphi_ref*costheta_ref*sinB_ref - cosA_ref*costheta_ref*sinphi_ref -
				sinA_ref*sintheta_ref + cosB_ref*(cosA_ref*sintheta_ref - costheta_ref*sinA_ref*sinphi_ref);
			// printf("%f, %ld \n", L_ref,L);
			#endif
			// L ranges from -sqrt(2) to +sqrt(2).  If it's < 0, the surface
			// is pointing away from us, so we won't bother trying to plot it.
			if (L > 0) {
				// test against the z-buffer.  larger 1/z means the pixel is
				// closer to the viewer than what's already plotted.
				if(ooz > zbuffer[yp][xp]) {
				zbuffer[yp][xp] = ooz;
				/* Here, the multiply by 8 has been transformed by removing some shift  */
				L >>= (SHIFT_QTY_U16-3);
				int luminance_index = L;

				// luminance_index is now in the range 0..11 (8*sqrt(2) = 11.3)
				// now we lookup the character corresponding to the
				// luminance and plot it in our output:
				output[yp][xp] = ".,-~:;=!*#$@"[luminance_index];
				}
			}
			}
		}

		/* Add new lines and terminates make it a simple printf */
		for (int j = 0; j < SCREEN_HEIGHT; j++) {
			// output[j][39] = ;
			output[j][39] = 10;
			if(j==(SCREEN_HEIGHT-1))
			{
				output[j][39] = 0;
			}
		}

		// now, dump output[] to the screen.
		// bring cursor to "home" location, in just about any currently-used
		// terminal emulation mode
		printf("\x1b[H");
		printf(output);


		A+=5;
		/* Check for rollovers to take care of array lookups */
		if(A > 314)
		{
			A -= 314;
		}
		B+=5;
		/* Check for rollovers to take care of array lookups */
		if(B > 314)
		{
			B -= 314;
		}
	}
	
}

void print_donut(void)
{	
	#ifdef DONUT_FLOATING
		donut_readable_float();
	#endif
	#ifdef DONUT_LOOKUP
		donut_readable_fixed_lookup();
	#endif
	#ifdef DONUT_FIXED_PT
		 donut_readable_all_fixed();
	#endif
	#ifdef DONUT_SHORT_FIXED
		donut_readable_all_fixed_u16();
	#endif

}



// void donut_original(void)
// {
// 	uint32_t k =0;
// 	float A=0,B=0,i,j,z[1760];
// 	char b[1760];
// 	printf("\x1b[2J");
// 	for(;;)
// 	{
// 		memset(b,32,1760);
// 		memset(z,0,1760*4);
// 		for(j=0;j<6.28;j+=0.30)
// 		{
// 			for(i=0;i<6.28;i+=0.10)
// 			{
// 				float c=sin(i),d=cos(j), f=sin(j), l=cos(i);
// 				float e= sin(A),g=cos(A), m=cos(B), n=sin(B);

// 				float h=d+2,D=1/(c* h*e+f*g+5),t=c*h*g-f*e;
// 				int x=40+30*D*(l*h*m-t*n),y=12+15*D*(l*h*n+t*m),o=x+80*y,N=8*((f*e-c*d*g)*m-c*d*e-f*g-l*d*n);
// 				if(22>y&&y>0&&x>0&&80>x&&D>z[o])
// 				{
// 					z[o]=D;
// 					b[o]=".,-~:;=!*#$@"[N>0?N:0];
// 				}
// 			}
// 		}
	
// 		printf("\x1b[H");
// 		for(k=0;k<1760+1;k++)
// 		{
// 			putchar(k%80?b[k]:10);
// 			A+=0.04;
// 			B+=0.02;
// 		}
// 	}
   
// }
#endif