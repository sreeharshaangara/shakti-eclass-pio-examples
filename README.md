# shakti-eclass-pio-examples
This repository holds some sample code for the Shakti e-class running on the Arty A35 board. Note that you will need to have the Shakti E-class bit-file programmed into the board for this to run. 

To build code using platform IO run the below, 

```
platformio run -e pinaka
```

To upload code, 

```
platformio run -e pinaka --target upload
```



### Special note: 

For running the donut example, please make sure your terminal window has atleast 40 height so the donut renders correctly. For TeraTerm, you can find this under Setup > Terminal > Terminal size