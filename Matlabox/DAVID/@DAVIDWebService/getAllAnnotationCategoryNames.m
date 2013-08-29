function xReturn = getAllAnnotationCategoryNames(obj)
%getAllAnnotationCategoryNames(obj)
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
    'getAllAnnotationCategoryNames', ...
    values,names,types,'document');
response = callSoapService( ...
    obj.endpoint, ...
    'urn:getAllAnnotationCategoryNames', ...
    soapMessage);
xReturn = parseSoapResponse(response);
