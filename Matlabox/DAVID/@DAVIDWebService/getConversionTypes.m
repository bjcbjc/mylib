function xReturn = getConversionTypes(obj)
%getConversionTypes(obj)
%
%     Input:
%   
%     Output:
%       return = (string)

% Build up the argument lists.
values = { ...
   };
names = { ...
   };
types = { ...
   };

% Create the message, make the call, and convert the response into a variable.
soapMessage = createSoapMessage( ...
    'http://service.session.sample', ...
    'getConversionTypes', ...
    values,names,types,'document');
response = callSoapService( ...
    obj.endpoint, ...
    'urn:getConversionTypes', ...
    soapMessage);
xReturn = parseSoapResponse(response);
