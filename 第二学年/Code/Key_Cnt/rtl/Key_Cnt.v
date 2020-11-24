module Key_Cnt(
    Clk,
    Rst_n,
    key1,
    key2,

    Rs232_Tx,
//    Tx_Done,
//    uart_state
);
		input Clk;
		input Rst_n;
		input key1;
		input key2;

		output wire Rs232_Tx;
//		output wire Tx_Done;
//		output wire uart_state;

		wire key_flag1, key1_state1;	
		wire key_flag2, key1_state2;	
		
		wire [7:0]Data_byte;
		wire Send_en;
	
   key_filter keyfilter1(
		.Clk			(Clk		),
		.Rst_n		(Rst_n		),
		.key_in		(key1		),
		.key_flag	(key_flag1	),
		.key_state	(key_state1	)
	);
	key_filter keyfilter2(
		.Clk			(Clk		),
		.Rst_n		(Rst_n		),
		.key_in		(key2		),
		.key_flag	(key_flag2	),
		.key_state	(key_state2	)
	);
	key_ctrl key_ctrl1(
		.Clk			(Clk		),
		.Rst_n		(Rst_n		),
		.key_flag0	(key_flag1	),
		.key_flag1	(key_flag2	),
		.key_state0	(key_state1	),
		.key_state1	(key_state2	),
		.Data_byte	(Data_byte),
		.Send_en	(Send_en)
	);
	
    uart_byte_tx uart_byte_tx1(
		.Clk			(Clk		),
		.Rst_n		(Rst_n		),
		.data_byte	(Data_byte	),
		.send_en	(Send_en),
		.baud_set	(3'd3),
		
		.Rs232_Tx	(Rs232_Tx	),
		.Tx_Done		(	),
		.uart_state	(	)
	);
	
endmodule