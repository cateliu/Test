module ctrl_uart(
	input 	Clk,
	input 	Rst_n,
//	input 	[31:0]T,
//	input 	wrt_en,
//	input 	Send_en;
//	output 	[7:0]data;
	output 	data_t
);
	wire Tx_Done;
	wire Send_en;
	wire[7:0] data;
	reg wrt_en;
	reg [31:0] cnt;
//	reg tx_sign;
	
//	always@(posedge Clk or negedge Rst_n)
//	if(!Rst_n)
//		tx_sign <= 1'b0;
//	else
//		tx_sign <= TX_Done;
		
	always@(posedge Clk or negedge Rst_n)begin
		if(!Rst_n)begin
			cnt 		<= 32'd0;
			wrt_en 	<= 1'b0;
		end
		else if(cnt == 32'h10_01_FE_FE)begin
			cnt		<= cnt +1'b1;
			wrt_en 	<= 1'b1;
		end
		else if(cnt == 32'h10_01_FE_FF)begin
			cnt		<= 32'd0;
			wrt_en	<= 1'b0;
		end
		else	begin
			cnt 		<= cnt +1'b1;
			wrt_en	<= 1'b0;
		end
	end
	ctrl ctrl1(
		.Clk			(Clk		),
		.Rst_n		(Rst_n	),
		.T				(cnt			),
		.TX_Done		(Tx_Done	),
		.wrt_en		(wrt_en	),
		.Send_en		(Send_en),
		.data			(data		)
	);
	
	uart_byte_tx uart_byte_tx1(
		.Clk			(Clk		),
		.Rst_n		(Rst_n	),
		.data_byte	(data		),
		.send_en		(Send_en	),
		.baud_set	(3'd4		),
		
		.Rs232_Tx	(data_t	),
		.Tx_Done		(Tx_Done	),
		.uart_state	()
	);
endmodule