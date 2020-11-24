module counter(Clk, Rst_n, Cin, T, c_en, M);
	
	input Clk;
	input Rst_n;
	input Cin;
	
	output reg[24:0] M;
	output reg c_en;
	output reg[31:0] T;
	
	reg [24:0] cnt;
	wire AB;
	reg A;
	reg B;
	// 连续两个输入量，如果第一个为0，第二个为1，那么输入信号为上升沿，AB=1；其余情况则为0；
	always@(posedge Clk or negedge Rst_n)
	if(!Rst_n) begin
		A <= 1'b0;
		B <= 1'b0;
	end
	else begin
		B <= Cin;
		A <= B;
	end
	
	assign AB = A&(!B);
	// 当上升沿到来时，重新计数，计数值cnt初始值为1.
	always@(posedge Clk or negedge Rst_n)
	if(!Rst_n)begin
		cnt <= 32'd1;
		c_en <= 1'b0;
		T <= 32'd0;
	end
	else if(AB === 1'b1)begin
		T <= cnt;
		cnt <= 32'd1;
		c_en <= 1'b1;
	end
	else if(cnt == 32'hFF_FF_FF_FF)begin
		cnt <= 32'd0;
	end
	else begin
		cnt <= cnt + 1'b1;
		c_en <= 1'b0;
		T	<= T;
	end
	
	// 计数方波信号的周期数
	always@(posedge Cin or negedge Rst_n)
	if(!Rst_n ) begin
		M <= 0;
	end
	else 
		M <= M + 1'b1;
	
	
	
endmodule